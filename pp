time ./raster_tuning_LIF -ng -v -inf - -n 80 20 -t 1e5 -seed 1 -mat net_100_0X65763652.txt -pr 1 -ps 0.008 -scee 0.005 -scei 0.005 -scie 0.009 -scii 0.009 --bin-save -o ./data/a.dat

real	1m26.416s
user	1m25.932s
sys	0m0.340s

icc xhp7
real	0m47.845s
user	0m47.060s
sys	0m0.160s


time bin/Release/vec_IFsimu --nE 80 --nI 20 --t 1e5 --net net_100_0X65763652.txt --pr 1 --ps 0.008 --scee 0.005 --scie 0.005 --scei 0.009 --scii 0.009 -o a.dat --isi-path=b.txt

real	0m22.122s
user	0m21.636s
sys	0m0.236s
xhp7
real	0m13.227s
user	0m13.084s
sys	0m0.124s

