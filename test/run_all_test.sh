#!/bin/bash
# copyright Scott Givan, The University of Missouri, July 6, 2012
#
#    This file is part of the RNA-seq Toolkit, or RST.
#
#    RST is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    RST is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with RST.  If not, see <http://www.gnu.org/licenses/>.
#
echo 'loading dependencies'

module load Python-shared
module load R-3.3.0-sharedlib
module load bowtie2-2.3.2
module load stringtie-1.3.0
module load HISAT2-2.1.0


echo 'Running all tests.'

for dir in `ls -1d TestData_*`
do
    echo ""
    echo "$dir"
    cd $dir
    ./setup_and_test.sh
    cd ..
done

