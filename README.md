# molecular-geometry-optimizer

# Compiling and Executing

- Compiling the program: `make all`

- Executing the program: `./moleculeGeometryOptimizer <input txt data file> <alpha electrons number> <beta electrons number> <SCF algorithm> <true/false optimize>`. For example: `./moleculeGeometryOptimizer ./data/C2H2.txt 5 5 DIIS flase > ./data/C2H2.txt.DIIS.out`

# Jmol commands for running animation

- Jmol data directory: `/private/var/folders/4b/cqxrc9r92jb1jc2qq33tn9jr0000gn/T/hsperfdata_xperrylinn/`

- load an xyz animation into Jmol: `load trajectory "animation.xyz"`

- copy command to .xyz animation files from data dir to Jmol dir: `cp data/*.xyz /private/var/folders/4b/cqxrc9r92jb1jc2qq33tn9jr0000gn/T/hsperfdata_xperrylinn/animations/`

- copy command to move from .xyz files from Jmol dir to data dir: `cp /private/var/folders/4b/cqxrc9r92jb1jc2qq33tn9jr0000gn/T/hsperfdata_xperrylinn/*.xyz ./data/Jmol_export_xyz/`

- animation sript (creates a sequence of jpg images):

 ```
 frame 1
 num_frames = getProperty("modelInfo.modelCount")
 for (var i = 1; i <= num_frames; i = i+1)
   var filename = "movie"+("00000"+i)[-4][0]+".jpg"
   write IMAGE 800 600 JPG @filename
   frame next
 end for
```

- write to .xyz file: `write my_molecule.xyz`

- cp .xyz from dir belonging to Jmol to data folder: `cp /private/var/folders/4b/cqxrc9r92jb1jc2qq33tn9jr0000gn/T/hsperfdata_xperrylinn/cyclohexane_chair.xyz ./data/cyclohexane_chair.xyz`