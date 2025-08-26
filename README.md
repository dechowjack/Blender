# Blender
Documentation and algorithm files for Blender [Dechow et al., 2025] to be submitted to Water Resources Research.
This is version 5.0 of the Blender algorithm, used to generate data for [Dechow et al., 2025] as well as my PhD Dissertation:

Merging remote sensing observations and land surface models to improve estimates of the spatial and temporal dynamics of snow water equivalent and surface density
http://rave.ohiolink.edu/etdc/view?acc_num=osu1723774509328696

Included are a sample job submission script [samplejob.txt] as well as the call script [RunTuolumne.sh] that runs the actual algorithm [Blender_Algorithm_v5.0.jl]
Unfortunately, Blender is written to run on a single pixel at a time. You will need to generate all input files on a per-pixel basis, and move them into named per-pixel directories
following the naming convention [Pix#]. A more detailed description of input parameters is available in my dissertation (as of August 26 2025) as well as [Dechow et al., 2025] # DOI will be added when published.
The actual algorithm is ran using a for/while loop to count over instances of calling [RunTuolumne.sh] with an int. input corresponding to the pixel you want to run. 
For example, to run Blender on the first pixel of your domain, create dir [domain/Pix1] and include input files [SWE.txt Precip.txt EnergyBalance.txt SCF.txt AirTemp.txt] inside
dir [Pix1], then in a bash environment run:

  sh RunTuolumne.sh 1

The algorithm script [Blender_Algorithm_v5.0.jl] will by default write all output files to input dir [Pix1]. Log files for the Julia IpOPT routine are generated but not saved unless a log file is defined during Julia call

see this line from [RunTuolumne.sh]:

  /apps/julia/1.1.1/bin/julia Blender_Algorithm_v5.0.jl . >log.txt

Further documentation may be added at a future date # -JLD 08/26/25
