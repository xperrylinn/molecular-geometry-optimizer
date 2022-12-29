# molecular-geometry-optimizer

# Compiling and Executing

- Compiling the program: `make all`

- Executing the program: `./moleculeGeometryOptimizer <input txt data file> <alpha electrons number> <beta electrons number>`. For example: `./moleculeGeometryOptimizer ./data/C2H2.txt 5 5 DIIS > ./data/C2H2.txt.DIIS.out`

# Jmol commands for running animation

- load trajectory "animation.xyz" 

- copy command to move to Jmol dir: `cp data/animation.xyz /private/var/folders/4b/cqxrc9r92jb1jc2qq33tn9jr0000gn/T/hsperfdata_xperrylinn/animation.xyz`

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