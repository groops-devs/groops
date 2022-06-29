# main
- New programs:     GriddedData2GriddedDataTimeSeries and GriddedDataTimeSeries2GriddedData.
- New programs:     MagneticField2GriddedData and Orbit2MagneticField.
- New class:        In MiscAccelerations: FromParametrization
- New option:       GnssAntennaDefinitionCreate: rename antennas.
- New option:       gnssReceiverGeneratorStationNetwork: inputfileClock.
- New option:       gnssReceiverGenerator: print preprocessing infos.
- New option:       GroupPrograms: silently and additional outputfileLog.
- Bugfix:           instrument files: empty files are now compatible to other instrument types.
- Bugfix:           gnssProcessingStep: uninitialized normalEquationInfo.
- Bugfix:           gnssProcessingStep: wrong counting of observations.
- Bugfix:           gnssProcessingStepForEachReceiverSeparately: variableReceiver was not set.
- Bugfix:           gnssProcessingStepResolveAmbiguities: for writing empty ambiguity file.
- Bugfix:           gnssParametrizationClocksModel: Fixed zero mean constraint.
- Bugfix:           gnssParametrizationLeoDynamicOrbits: in parallel excecution.
- Bugfix:           gnssParametrizationKinematicPositions: in parallel excecution.
- Bugfix:           gnssTransmitter: noAntennaPatternFound->ignoreObservation not working correctly.
- Bugfix:           gnssReceiver: Simulating GLONASS ambiguities now correctly considers frequency channel.
- Bugfix:           sp3Format2Orbit: no/invalid orbit positions/velocities are now excluded.
- Bugfix:           Conversion of GRACE L1B/L1A data: revised source code.
- Bugfix:           loopFileAscii: Fixed uninitialized variable that could lead to the loop ending prematurely.
- Bugfix:           GnssAntex2AntennaDefinition: Fixed handling of frequency RMS blocks.
- Other:            File GriddedDataTimeSeries: includes now the last epoch; interval [...] instead of [...).
- Other:            File TimeSplinesGravityfield: includes now the last epoch; interval [...] instead of [...).
- Other:            Removed inputfileGlobal option.
- Other:            GnssAttitude2Orbex: can now handle different sampling per satellite.
- Other:            GnssRinexNavigation2OrbitClock/RinexObservation2GnssReceiver: Added basic support for RINEX v4.00.
- Other:            gnssParametrization*DynamicOrbits: integration starts and ends with first/last valid epoch.
- Other:            GnssLowEarthOrbiter: createTracks() before removing outlier epochs leads to less track splits.

# Release 2021-09-06
- Interface change: Complete redesign of GnssProcessing to make usage a little bit easier and more flexible.
    - Direct use of orbits without integrating variational equations in case of fixed transmitters (e.g., PPP).
    - New class to add flexible parametrizations to the normal equation system.
    - New class to select transmitters/receivers for each parametrization.
    - Unified all transmitter classes into single class and merged all transmitter data and metadata into one folder at https://ftp.tugraz.at/outgoing/ITSG/groops/data/gnss/.
    - Example scenarios with config files at https://ftp.tugraz.at/outgoing/ITSG/groops/scenario/.
    - Updated and expanded documentation and cookbooks to reflect all GNSS-related changes.
- New program:      InstrumentAccelerometer2ThermosphericDensity: Estimate neutral density from accelerometer data.
- New class:        In Thermosphere: new model nrlmsis2
- New class:        In Condition: Matrix to evaluate matrix elements.
- New class:        In PlotMapProjection: added Mollweide map projection.
- Bugfix:           FileSatelliteModel: removed if-statement for shaded plates, not necessary when applying the algorithm following Sentman 1961
- Other:            Expression parser: constants are now defined with brackets, e.g pi().
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
