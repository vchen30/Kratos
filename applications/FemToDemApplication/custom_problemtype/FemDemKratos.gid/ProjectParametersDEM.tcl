
proc WriteProjectParametersDEM { basename dir problemtypedir } {


    set filename [file join $dir ProjectParametersDEM.json]
    set FileVar  [open $filename w]

    puts $FileVar "{"
    puts $FileVar "   \"Dimension\"                      : [GiD_AccessValue get gendata Dimension],"
    puts $FileVar "   \"PeriodicDomainOption\"           : false,"
    puts $FileVar "   \"BoundingBoxOption\"              : false,"
    puts $FileVar "   \"AutomaticBoundingBoxOption\"     : false,"
    puts $FileVar "   \"BoundingBoxEnlargementFactor\"   : 1.0,"
    puts $FileVar "   \"BoundingBoxStartTime\"           : 0.0,"
    puts $FileVar "   \"BoundingBoxStopTime\"            : 1000.0,"
    puts $FileVar "   \"BoundingBoxMaxX\"                : 10,"
    puts $FileVar "   \"BoundingBoxMaxY\"                : 10,"
    puts $FileVar "   \"BoundingBoxMaxZ\"                : 10,"
    puts $FileVar "   \"BoundingBoxMinX\"                : -10,"
    puts $FileVar "   \"BoundingBoxMinY\"                : -10,"
    puts $FileVar "   \"BoundingBoxMinZ\"                : -10,"
    puts $FileVar ""
    puts $FileVar "   \"dem_inlet_option\"               : false,"
    puts $FileVar "   \"GravityX\"                       : [GiD_AccessValue get gendata GravityX],"
    puts $FileVar "   \"GravityY\"                       : [GiD_AccessValue get gendata GravityY],"
    puts $FileVar "   \"GravityZ\"                       : [GiD_AccessValue get gendata GravityZ],"
    puts $FileVar ""
    puts $FileVar "   \"EnergyCalculationOption\"        : false,"
    puts $FileVar "   \"VelocityTrapOption\"             : false,"
    puts $FileVar "   \"RotationOption\"                 : true,"
    puts $FileVar "   \"CleanIndentationsOption\"        : false,"
    puts $FileVar "   \"RemoveBallsInEmbeddedOption\"    : true,"
    puts $FileVar ""
    puts $FileVar "   \"DeltaOption\"                    : \"Absolute\","
    puts $FileVar "   \"SearchTolerance\"                : 0.0,"
    puts $FileVar "   \"AmplifiedSearchRadiusExtension\" : 0.0,"
    puts $FileVar "   \"ModelDataInfo\"                  : false,"
    puts $FileVar "   \"VirtualMassCoefficient\"         : 1.0,"
    puts $FileVar "   \"RollingFrictionOption\"          : false,"
    puts $FileVar "   \"ContactMeshOption\"              : false,"
    puts $FileVar "   \"OutputFileType\"                 : \"Binary\","
    puts $FileVar "   \"Multifile\"                      : \"multiple_files\","
    puts $FileVar "   \"ElementType\"                    : \"SphericPartDEMElement3D\","
    puts $FileVar ""
    puts $FileVar "   \"IntegrationScheme\"              : \"Symplectic_Euler\","
    puts $FileVar "   \"AutomaticTimestep\"              : false,"
    puts $FileVar "   \"DeltaTimeSafetyFactor\"          : 1.0,"
    puts $FileVar "   \"MaxTimeStep\"                    : [GiD_AccessValue get gendata MaxTimeStep],"
    puts $FileVar "   \"FinalTime\"                      : 1.0,"
    puts $FileVar "   \"ControlTime\"                    : 4.0,"
    puts $FileVar "   \"NeighbourSearchFrequency\"       : 1,"




    puts $FileVar "}"
    close $FileVar
}