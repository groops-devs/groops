![GROOPS Logo](https://github.com/groops-devs/groops/blob/main/docs/html/static/groops_banner.png)

The Gravity Recovery Object Oriented Programming System (GROOPS) is a software toolkit written in C++
that enables the user to perform core geodetic tasks.
Key features of the software include gravity field recovery from satellite and terrestrial data,
the determination of satellite orbits from global navigation satellite system (GNSS) measurements,
and the processing of GNSS constellations and ground station networks.

Most tasks and algorithms are (optionally) parallelized through the Message Passing Interface (MPI), thus
the software enables a smooth transition from single-CPU desktop computers to large distributed
computing environments for resource intensive tasks.

For an easy and intuitive setup of complex workflows, GROOPS contains a graphical
user interface where configuration files can be created and edited.

- [Citing GROOPS](#citing-groops)
- [Installation](#installation)
- [Getting Started](#getting-started)
- [Contributing](#contributing)
- [License](#license)
- [Contributors](#contributors)

## Citing GROOPS

If you use data sets computed with GROOPS in a publication or publish the data itself,
please cite our [reference paper](https://doi.org/10.1016/j.cageo.2021.104864):

*Mayer-Guerr, T., Behzadpour, S., Eicker, A., Ellmer, M., Koch, B., Krauss, S., Pock, C., Rieser, D., Strasser, S., Suesser-Rechberger, B., Zehentner, N.,  Kvas, A. (2021). GROOPS: A software toolkit for gravity field recovery and GNSS processing. Computers & Geosciences, 104864. https://doi.org/10.1016/j.cageo.2021.104864*

```
@article{Mayer-Gurr2021,
  author = {Mayer-Guerr, Torsten and Behzadpour, Saniya and Eicker, Annette and Ellmer, Matthias and Koch, Beate and Krauss, Sandro and Pock, Christian and Rieser, Daniel and Strasser, Sebastian and Suesser-Rechberger, Barbara and Zehentner, Norbert and Kvas, Andreas},
  doi = {https://doi.org/10.1016/j.cageo.2021.104864},
  issn = {0098-3004},
  journal = {Computers & Geosciences},
  keywords = {GNSS processing,Gravity field recovery,Orbit determination},
  pages = {104864},
  title = {{GROOPS: A software toolkit for gravity field recovery and GNSS processing}},
  url = {https://www.sciencedirect.com/science/article/pii/S009830042100159X},
  year = {2021}
}
```

## Installation

GROOPS is written in C++ and contains some legacy Fortran code.
To enable an intuitive interaction with the software, GROOPS includes a
graphical user interface (GUI).
The GUI is also written in C++ and depends on the Qt toolkit.

A detailed installation guide for Microsoft Windows and various Linux distributions can be found
on the [Installation page](https://github.com/groops-devs/groops/blob/main/INSTALL.md).

## Getting Started

After a successful installation our [Documentation](https://groops-devs.github.io/groops/html/index.html)
is the perfect way to get familiar with the different features of GROOPS.

GROOPS depends on data files such as Earth rotation, Love numbers, and wavelet coefficients.
An initial data set that is regularly updated is available on [our FTP server](https://ftp.tugraz.at/outgoing/ITSG/groops/).
You can choose between downloading the data directory or
a single [zip file](https://ftp.tugraz.at/outgoing/ITSG/groops/data.zip) with the same content.

## Contributing

We appreciate all contributions such as improving the documentation, reporting or fixing bugs,
implementing new features. Answering user questions in the
[Discussions](https://github.com/groops-devs/groops/discussions) section is another great way
of contributing to the GROOPS community.

If you encounter a bug, please let us know by [filing an issue](https://github.com/groops-devs/groops/issues).
Please include as much information as possible on how to reproduce the bug
and about your software environment (operating system, compiler version, GROOPS version).

If you want to provide a bug fix or implement a new features,
please get in contact with us in the [Discussions](https://github.com/groops-devs/groops/discussions)
before you start coding.

## License

GROOPS is licensed under GPLv3, as found in the [LICENSE](https://github.com/groops-devs/groops/blob/main/LICENSE) file.
This license applies to all files in the repository unless otherwise indicated.

Information about external source code contained in the repository which is licensed differently can be found in the
[corresponding README](https://github.com/groops-devs/groops/blob/main/source/external/README.md).

## Contributors

Parts of GROOPS originate from developments in the Astronomical, Physical and Mathematical Geodesy Group
at the University of Bonn, Germany.
Since 2010 it is developed and maintained at Graz University of Technology, Austria.

Here is a list of current and past contributors:

[Torsten Mayer-Guerr](https://github.com/tmayerguerr), Annette Eicker, Daniel Rieser, Norbert Zehentner,
Christian Pock, [Matthias Ellmer](https://github.com/x49), Beate Koch, [Andreas Kvas](https://github.com/akvas), Saniya Behzadpour,
[Sebastian Strasser](https://github.com/sestras), Sandro Krauss, Barbara Suesser-Rechberger,
[Patrick Dumitraschkewitz](https://github.com/zhedumi), Felix Oehlinger
