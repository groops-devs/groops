# main
- New program:      InstrumentAccelerometer2ThermosphericDensity: Estimate neutral density from accelerometer data.
- New class:        In Thermosphere: new model nrlmsis2
- New class:        In Condition: Matrix to evaluate matrix elements.
- New class:        GnssParametrizationTransmitterQzss: QZSS constellation support in GNSS processing.
- New class:        In PlotMapProjection: added Mollweide map projection.
- Bugfix:           FileSatelliteModel: removed if-statement for shaded plates, not necessary when applying the algorithm following Sentman 1961
- Other:            Additional constants in the mathematical parser like speed of light c().

# Release 2021-02-02
- Interface change: GnssProcessing, GnssSimulateReceiver: Removed intervals (use program within LoopPrograms instead).
- Interface change: SimulateStarCameraGnss: Full reimplementation with interface change.
                    Added support for all known attitude modes used by GPS, GLONASS, Galileo, BeiDou, and QZSS. Now requires GnssAttitudeInfo file.
- Interface change: Renamed program KalmanStaticTemporalNormals to NormalsBuildShortTimeStaticLongTime.
- New program:      GnssAttitudeInfoCreate: Creates attitude info file used by SimulateStarCameraGnss.
- New program:      PreprocessingDualSst: Analyze GRACE-FO KBR and LRI together.
- New class:        In Observation: DualSstVariational to use GRACE-FO KBR and LRI together.
- New class:        In ParametrizationGravity: LinearTransformation: Gravity field parametrization based on the linear transformation of another parametrizationGravity.
- New option:       LoopPrograms: processCountPerIteration (when running the loop on multiple processes), parallelLog (output to screen/log files from all processes).
- New option:       IfPrograms: elsePrograms (executed if condition evaluates to false).
- New option:       GroupPrograms: catchErrors (prevents program termination on error and optionally runs additional programs, i.e. try-catch).
- Bugfix:           Orbit2Kepler: Fixed angular output values (DEG2RAD -> RAD2DEG).
- Bugfix:           GnssClockRinex2InstrumentClock: 9-character identifier field width is now used starting from v3.04, not (incorrectly) from v3.00.
- Bugfix:           SphericalHarmonicsFilterMatrix: Input coefficient vector is now sorted correctly into filter matrix numbering.
- Bugfix:           MatrixDistributed: choleskyInverse(): Fixed a bug with sparse matrices.
- Bugfix:           Rectangular grids with one row or column (i.e. parallels or meridians) are now handled correctly.
- Bugfix:           InstrumentEstimateEmpiricalCovariance: Computation of autocovariance now works as expected.
- Bugfix:           Parallel: Multiple bugfixes and improvements for better support of different MPI implementations.
- Other:            Gnss: Updated BeiDou signal definition according to RINEX 3.05 and added support for BeiDou composite types.
- Other:            Sp3Format2Orbit: Added support for SP3d format.
- Other:            LoopPrograms: continueAfterError now works in parallel execution.
- Other:            Improved CMake installation process (see updated INSTALL.md). Now supports parallel compilation and install target.

# Release 2020-11-12
- Initial release
