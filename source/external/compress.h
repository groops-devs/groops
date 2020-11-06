/********************************************************************************************
* Public domain implementation of the LZW algorithm. See https://github.com/vapier/ncompress.
*/

#ifndef GROOPS_UNCOMPRESS_H
#define GROOPS_UNCOMPRESS_H


#include <string.h>
#include <sstream>
#include <fstream>

#define IBUFSIZ 2048 /* Defailt input buffer size							*/
#define OBUFSIZ 2048 /* Default output buffer size							*/

#define MAGIC_1 (unsigned char)'\037' /* First byte of compressed file				*/
#define MAGIC_2 (unsigned char)'\235' /* Second byte of compressed file				*/
#define BIT_MASK 0x1f                 /* Mask for 'number of compresssion bits'		*/
                                      /* Masks 0x20 and 0x40 are free.  				*/
                                      /* I think 0x20 should mean that there is		*/
                                      /* a fourth header byte (for expansion).    	*/
#define BLOCK_MODE 0x80               /* Block compresssion if table is full and		*/
                                      /* compression rate is dropping flush tables	*/
                                      /* the next two codes should not be changed lightly, as they must not	*/
                                      /* lie within the contiguous general code space.						*/
#define FIRST 257                     /* first free entry 							*/
#define CLEAR 256                     /* table clear output code 						*/

#define INIT_BITS 9 /* initial number of bits/code */

#define HBITS 17   /* 50% occupancy */
#define HSIZE (1 << HBITS)

#define BITS 16

#define MAXCODE(n) (1L << (n))

union bytes {
  long word;
  struct
  {
    unsigned char b4;
    unsigned char b3;
    unsigned char b2;
    unsigned char b1;
  } bytes;
};

#define input(b, o, c, n, m)                       \
  {                                                \
    unsigned char *p = &(b)[(o) >> 3];             \
    (c) = ((((long)(p[0])) | ((long)(p[1]) << 8) | \
            ((long)(p[2]) << 16)) >>               \
           ((o)&0x7)) &                            \
          (m);                                     \
    (o) += (n);                                    \
  }

/*
 * Decompress stdin to stdout.  This routine adapts to the codes in the
 * file building the "string" table on-the-fly; requiring no table to
 * be stored in the compressed file.  The tables used herein are shared
 * with those of the compress() routine.  See the definitions above.
 */

void decompress(std::ifstream& infile, std::stringstream &ss)
{
  unsigned char inbuf[IBUFSIZ + 64];    /* Input buffer									*/
  unsigned char outbuf[OBUFSIZ + 2048]; /* Output buffer								*/

  long int htab[HSIZE];
  unsigned short codetab[HSIZE];

  int insize = 0;
  int rsize;
  while (insize < 3 && (rsize = infile.readsome((char*)&inbuf + insize, IBUFSIZ)) > 0)
    insize += rsize;

  if(insize < 3 || inbuf[0] != MAGIC_1 || inbuf[1] != MAGIC_2)
  {
    if(rsize < 0)
      throw(std::ios_base::failure("uncompress: unable to read input file"));

    if (insize > 0)
      throw(std::ios_base::failure("uncompress: input file does not seem to be in compressed format."));
  }

  int maxbits = inbuf[2] & BIT_MASK;
  int block_mode = inbuf[2] & BLOCK_MODE;
  long int maxmaxcode = MAXCODE(maxbits);

  if (maxbits > BITS)
    throw(std::ios_base::failure("uncompress: number of compression bits in input file does not match algorithm."));

  int n_bits;
  long int maxcode = MAXCODE(n_bits = INIT_BITS) - 1;
  int bitmask = (1 << n_bits) - 1;
  long int oldcode = -1;
  int finchar = 0;
  int outpos = 0;
  int posbits = 3 << 3;

  long int free_ent = ((block_mode) ? FIRST : 256);

  memset(codetab, 0, 256); /* As above, initialize the first 256 entries in the table. */
  long int code;
  for (code = 255; code >= 0; --code)
    ((unsigned char *)(htab))[code] = (unsigned char)code;

  do
  {
  resetbuf:;
    {
      int o = posbits >> 3;
      int e = o <= insize ? insize - o : 0;

      for (int i = 0; i < e; ++i)
        inbuf[i] = inbuf[i + o];

      insize = e;
      posbits = 0;
    }

    if (insize < (int)sizeof(inbuf) - IBUFSIZ)
    {
      if ((rsize = infile.readsome((char*)&inbuf + insize, IBUFSIZ)) < 0)
        throw(std::ios_base::failure("uncompress: unable to read input file"));

      insize += rsize;
    }

    int inbits = ((rsize > 0) ? (insize - insize % n_bits) << 3 : (insize << 3) - (n_bits - 1));

    while (inbits > posbits)
    {
      if (free_ent > maxcode)
      {
        posbits = ((posbits - 1) + ((n_bits << 3) - (posbits - 1 + (n_bits << 3)) % (n_bits << 3)));

        ++n_bits;
        if (n_bits == maxbits)
          maxcode = maxmaxcode;
        else
          maxcode = MAXCODE(n_bits) - 1;

        bitmask = (1 << n_bits) - 1;
        goto resetbuf;
      }

      input(inbuf, posbits, code, n_bits, bitmask);

      if (oldcode == -1)
      {
        if (code >= 256)
          throw(std::ios_base::failure("uncompress: input seems to be corrupt."));

        outbuf[outpos++] = (unsigned char)(finchar = (int)(oldcode = code));
        continue;
      }

      if (code == CLEAR && block_mode)
      {
        memset(codetab, 0, 256); //clear_tab_prefixof();
        free_ent = FIRST - 1;
        posbits = ((posbits - 1) + ((n_bits << 3) -
                                    (posbits - 1 + (n_bits << 3)) % (n_bits << 3)));
        maxcode = MAXCODE(n_bits = INIT_BITS) - 1;
        bitmask = (1 << n_bits) - 1;
        goto resetbuf;
      }

      long int incode = code;
      unsigned char *stackp = ((unsigned char *)&(htab[HSIZE - 1]));

      if (code >= free_ent) /* Special case for KwKwK string.	*/
      {
        if (code > free_ent)
        {
          unsigned char *p;

          posbits -= n_bits;
          p = &inbuf[posbits >> 3];

          fprintf(stderr, "insize:%d posbits:%d inbuf:%02X %02X %02X %02X %02X (%d)\n", insize, posbits,
                  p[-1], p[0], p[1], p[2], p[3], (posbits & 07));
          fprintf(stderr, "uncompress: corrupt input\n");
          throw(std::exception()); //exit(1); // read error
        }

        *--stackp = (unsigned char)finchar;
        code = oldcode;
      }

      while ((long int)code >= (long int)256)
      { /* Generate output characters in reverse order */
        *--stackp = ((unsigned char *)(htab))[code];
        code = codetab[code];
      }
      *--stackp = (unsigned char)(finchar = ((unsigned char *)(htab))[code]);

      /* And put them out in forward order */
      {
        int i;

        if (outpos + (i = (((unsigned char *)&(htab[HSIZE - 1])) - stackp)) >= OBUFSIZ)
        {
          do
          {
            if (i > OBUFSIZ - outpos)
              i = OBUFSIZ - outpos;

            if (i > 0)
            {
              memcpy(outbuf + outpos, stackp, i);
              outpos += i;
            }

            if (outpos >= OBUFSIZ)
            {
              for (int k = 0; k < outpos; k++)
                ss << outbuf[k];

              outpos = 0;
            }
            stackp += i;
          } while ((i = (((unsigned char *)&(htab[HSIZE - 1])) - stackp)) > 0);
        }
        else
        {
          memcpy(outbuf + outpos, stackp, i);
          outpos += i;
        }
      }

      if ((code = free_ent) < maxmaxcode) /* Generate the new entry. */
      {
        codetab[code] = (unsigned short)oldcode;
        ((unsigned char *)(htab))[code] = (unsigned char)finchar;
        free_ent = code + 1;
      }
      oldcode = incode; /* Remember previous code.	*/
    }
  } while (rsize > 0);

  for (int k = 0; k < outpos; k++)
    ss << outbuf[k];
}

/*****************************************************************
 * Algorithm:
 *   Modified Lempel-Ziv method (LZW).  Basically finds common
 *   substrings and replaces them with a variable size code.  This is
 *   deterministic, and can be done on the fly.  Thus, the decompression
 *   procedure needs no input table, but tracks the way the table was built.
 */
std::string decompress_file(const char *fileName)
{
  std::ifstream infile(fileName, std::ios::in|std::ios::binary);
  if(!infile.good())
    throw(std::runtime_error("error by opening file"));

  std::stringstream ss;
  decompress(infile, ss);

  return ss.str();
}

#endif
