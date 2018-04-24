#!/bin/bash
#
#   Copyright Scott A. Givan, University of Missouri, July 6, 2012.
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
./reset_test
export HISAT_VERSION_BEING_TESTED=2.1.0
TESTDATA_TEST_LOG="test_${HISAT_VERSION_BEING_TESTED}.log";
t/versionless_setup.sh |& tee $TESTDATA_TEST_LOG
t/run_test.sh          |& tee --append $TESTDATA_TEST_LOG
t/test.t $HISAT_VERSION_BEING_TESTED $@ |& tee --append $TESTDATA_TEST_LOG
