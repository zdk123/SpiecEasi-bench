#!/bin/bash
# Build shim dylibs that stand in for R's libRlapack (and OpenBLAS) so ./rbench can
# swap LAPACK/BLAS backends via DYLD_LIBRARY_PATH without modifying the R install.
set -e
cd "$(dirname "$0")"
rlapack=$(otool -D "$(Rscript -e 'cat(R.home("lib"))')/libRlapack.dylib" | tail -1)
openblas=/opt/homebrew/opt/openblas/lib/libopenblas.0.dylib
mkdir -p openblas accelerate
clang -dynamiclib shim.c -o openblas/libRlapack.dylib -install_name "$rlapack" -Wl,-reexport_library,$openblas
clang -dynamiclib shim.c -o accelerate/libRlapack.dylib -install_name "$rlapack" -Wl,-reexport_framework,Accelerate
clang -dynamiclib shim.c -o accelerate/libopenblas.0.dylib -install_name $openblas -Wl,-reexport_framework,Accelerate
