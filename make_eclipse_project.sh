#!/bin/sh
### Set BUILD_DIR and SOURCE_DIR appropriately
set -e
cd ~/Workspace/eclipse
rm -rf bsp-eclipse
mkdir bsp-eclipse
cd bsp-eclipse
cmake -G"Eclipse CDT4 - Unix Makefiles" ~/Workspace/bsp/trajopt
