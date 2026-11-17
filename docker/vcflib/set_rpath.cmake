# cmake -DF=<elf file> -P set_rpath.cmake: give a rebuilt vcflib file the RPATH the bioconda
# package's files have (docker/vcflib/build_vcflib.sh).
file(RPATH_SET FILE "${F}" NEW_RPATH "$ORIGIN/../lib")
