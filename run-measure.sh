#! /bin/sh

memtime="test-example7-mhap-lambda.sh.memtime"
/usr/bin/time --format "Command line: %C\nReal time: %e s\nCPU time: -1.0 s\nUser time: %U s\nSystem time: %S s\nMaximum RSS: %M kB\nExit status: %x" --quiet -o $memtime \
	./example7-mhap-lambda.sh

memtime="test-example8-mhap-ecoli_map006.sh.memtime"
/usr/bin/time --format "Command line: %C\nReal time: %e s\nCPU time: -1.0 s\nUser time: %U s\nSystem time: %S s\nMaximum RSS: %M kB\nExit status: %x" --quiet -o $memtime \
	./example8-mhap-ecoli_map006.sh

memtime="test-example9-mhap-ecoli_pacbio_160x.sh.memtime"
/usr/bin/time --format "Command line: %C\nReal time: %e s\nCPU time: -1.0 s\nUser time: %U s\nSystem time: %S s\nMaximum RSS: %M kB\nExit status: %x" --quiet -o $memtime \
	./example9-mhap-ecoli_pacbio_160x.sh
