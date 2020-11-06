## External source code used in GROOPS

While GROOPS is intended to be a standalone software package, some functionality depends on external software.
We list external source code contained in the GROOPS repository with references and licenses here.

### Horizontal Wind Model 2014 (HWM14)

* Reference: Drob, D. P., Emmert, J. T., Meriwether, J. W., Makela, J. J., Doornbos, E., Conde, M., Hernandez, G., Noto, J., Zawdie, K. A., McDonald, S. E., et al. (2015), An update to the Horizontal Wind Model (HWM): The quiet time thermosphere, Earth and Space Science, 2, 301– 319, https://doi.org/10.1002/2014EA000089
* Files: `hwm/hwm14.f90`, `hwm/README.txt`
* License: [CC BY-NC-ND 4.0](https://creativecommons.org/licenses/by-nc-nd/4.0)

### The Jacchia-Bowman 2008 Empirical Thermospheric Density Model (JB2008)

* Reference: Bowman, B. R., Tobiska, W. K., Marcos, F. A., Huang, C. Y., Lin, C. S., & Burke, W. J. (2008). A new empirical thermospheric density model JB2008 using new solar and geomagnetic indices. In AIAA/AAS Astrodynamics Specialist Conference and Exhibit. https://doi.org/10.2514/6.2008-6438
* Files: `jb2008/JB2008.f`, `jb2008/README.txt`
* License: Space Environment Technologies License and Warranty Agreement (http://sol.spacenvironment.net/jb2008/License.html)

### International Geomagnetic Reference Field (IGRF)

* Reference: Thébault, E., Finlay, C.C., Beggan, C.D. et al. International Geomagnetic Reference Field:
  the 12th generation. Earth Planet Sp 67, 79 (2015). https://doi.org/10.1186/s40623-015-0228-9
* Files: `igrf/igrf13.f`
* License: none (included with permission of the authors)

### International Earth Rotation and Reference Systems Service (IERS) Conventions software collection

* Reference: IERS Conventions (2010). Gérard Petit and Brian Luzum (eds.). (IERS Technical Note ; 36)
  Frankfurt am Main: Verlag des Bundesamts für Kartographie und Geodäsie, 2010. 179 pp., ISBN 3-89888-989-6
* Files: `iers/*` except `iers.h`
* License: IERS Conventions Software License
