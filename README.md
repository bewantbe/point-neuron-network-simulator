Accurate Point Neuronal Network Simulator
=========================================
(Note: **NOT** for artificial neural network)

This program aims to provide an accurate simulator and a portable code framework for point spiking neuron model networks.

A handy GNU Octave/Matlab interface to the core simulator is also provided.

The accuracy is achieved by breaking the simulation time steps (dt) according to the (e.g. spike) event timings. So that the dynamics inside the sub-time-steps are smooth and the classical ODE solves (notably RK4) work well.

Currently neuron models are (parenthesis: Name used inside the code):

  * Leaky Integrate-and-Fire models with exponential decay conductance (LIF-G).
  * Leaky Integrate-and-Fire models with "alpha function" conductance (LIF-GH).
  * Hodgkin Huxley neuron model with exponential decay conductance (HH-G).
  * Hodgkin Huxley neuron model with "alpha function" conductance (HH-GH).
  * Hodgkin Huxley neuron model with continous (but mostly around spike) synaptic interaction (HH-GH-cont-syn).
  * Jump current neuron model (IF-jump).
  * Hawkes model (Hawkes-GH).

Optional modification to the neuron models:

  * Add an alternating current to these models.(--current-sine-amp, --current-sine-freq)
  * Add an constant current to these models.(--extI)
  * Add synaptic delay. The delay must larger than dt.(--synaptic-delay or --synaptic-net-delay)
  * For HH type neuron, possible to use the peak or falling-threshold-passing as the spike time. (e.g. HH-PT-GH, HH-PT-G)
  * Note: some combinations of neuron model and modification are not yet implemented.

See [doc/neuron_models.pdf](doc/neuron_models.pdf) for model details of neuron models.

Features:

  * You can provide external input events in a file or just let the program generate (Poisson input) for you.
  * Input events can be excititory or/and inhibitory, and the strength can be specified at each timing. Each line of the input event file is "id time strength", the "strength" term is optional (default to --ps), negative strength means it is inhibitory.
  * Network specification (--net) supports adjacency matrix in either full form or sparse form (each line is "i j a\_{ij}"), and the final coupling strength is a\_{ij} times one of scee scie scei scii.
  * Neuron indexes are 1-based, but inside the code it is 0-based.

Refer to `bin/gen_neu --help` for full command line option help (after compilation).

### Note:

   Code still in developing stage, use with care.


Comparison to other simulators
------------------------------

* [NEST](http://www.nest-simulator.org/)

    NEST is a full featured point neuron simulator. While this one (APNNS) has much less neuron models and tunable settings. The function of APNNS is extended by directly modifying the source code.

    NEST is a big project. APNNS is rather small, roughly about 5 thousand lines of C++ code.
    
    NEST use a high order ODE solver for each neurons (by default it is Runge–Kutta–Fehlberg method (aka rkf45 or ode45()) through GNU Scientific Library). But the interactions between neurons has only order one accuracy, because the interactions is performed at the boundary of time step. APNNS, in contrast, will deal with these spike timings accurately.

* [NEURON](https://www.neuron.yale.edu/neuron/)

    NEURON is a more detailed/physiology simulator. APNNS only simulate "point" neurons.

* [Brian](http://briansimulator.org/)

    A clock-driven simulator, includes several integrators both for ODE and SDE. The spikes can only occur on the time grid hence the overall accuracy is O(dt).

    The dynamical equations are input as code string and then be translate to native code (e.g. by cpython) for computation. This makes the simulator very flexible and extensible while been fast.


Build from source
-----------------

### Build under Linux

You will need to install libraries Eigen and Boost.ProgramOptions first.
Then simply `make`, you will get an executable `bin/gen_neu`.


### Build the Matlab interface<a name="build-matlab"></a>

There are matlab scripts in `mfile/`, notably the interface `mfile/gen_neu.m`.

To use the matlab interface, you need to compile `mfile/BKDRHash.c` and `mfile/randMT19937.cpp`.

In matlab, use `mfile/` as your working directory.

  * compile `mfile/BKDRHash.c`

		mex -largeArrayDims BKDRHash.c
	
  * compile `mfile/randMT19937.cpp`

	Using GCC, under Linux

		mex CXXFLAGS="\$CXXFLAGS -std=c++11" CXXOPTIMFLAGS="-O3 -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" randMT19937.cpp

	Using MSVC, under Windows

		mex COMPFLAGS="$COMPFLAGS" OPTIMFLAGS="$OPTIMFLAGS /fopenmp" LINKFLAGS="$LINKFLAGS /fopenmp" randMT19937.cpp
	
	Using GCC, Octave, under Linux. (In command line)

		CXXFLAGS="-O3 -fopenmp -std=c++11"  LDFLAGS=" -fopenmp" mkoctfile --mex randMT19937.cpp


### 在 Windows 下用 MSVC 编译 (Build under Windows using MSVC)

示例环境：

  Win 10, Microsoft Visual Studio Community 2015

步骤：

1. 下载 boost

		http://www.boost.org/users/download/
		点 Prebuilt windows binaries.
		(例如 https://sourceforge.net/projects/boost/files/boost-binaries/1.60.0/)
		VS -> Help -> About 查看 VS 的版本。 比如 2015 是 version 14.0
		下载(32bit 通用性好一些) boost_1_60_0-msvc-14.0-32.exe
		安装，按默认路径。

2. 下载 Eigen

		http://eigen.tuxfamily.org/
		选择 "latest stable release" 
		（例如 http://bitbucket.org/eigen/eigen/get/3.2.8.zip）
    
		解压复制到编译系统文件夹
			C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\include
		(如果找不到 VC 的 inclulde 目录，那么可以打开 common_header.h, 右键查看 #include <vector>, 找到所在的路径。)
		重命名文件夹为 eigen3
		(也就是应该有文件 C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\include\eigen3\Eigen\Dense)
	
		也可以放到任意放一个地方（路径最好不要含中文）, 然后添加 Include Directories

3. 用 MSVC 打开 point-neuron-network-simulator.sln
  或手动新建工程，添加需要的文件。
  
  * 设置 boost

			手动在工程中添加 include 目录（若用默认路径安装，可能不需要设定）：
			在工程属性 (菜单项目->属性) 里设定 (C/C++ -> 常规 -> 附加包含目录，链接器 -> 常规 -> 附加库目录)
				Include Directories 为 C:\local\boost_1_60_0
				Library Directories 为 C:\local\boost_1_60_0\lib32-msvc-14.0
			如果是 64bit 工程，则
				Library Directories 为 C:\local\boost_1_60_0\lib64-msvc-14.0

  * 设置 Eigen

		如库文件没有放到编译系统文件夹，请自行设定 Include Directories。

	测试用命令行参数

		--neuron-model HH-GH --net - --dt 0.03125 --t 1e3 --pr 1 --ps 0.05  --seed 1

	没报错就是编译正常了。注意要用 Release 模式编译以保证性能，以及选择 x86 ()。
	
	把生成的可执行程序(通常在 Release/ 目录下) 重命名成 gen_neu.exe 并放到 mfile/ 目录下。

Usage
-----

### Using the core simulator `bin/gen_neu`

(Note: compile the program first)

Enter `bin/gen_neu --help` to see command line help.

Usage example:

		bin/gen_neu --neuron-model HH-GH --net - --nE 1 --nI 0 --t 1e3 --dt 0.03125 --pr 1 --ps 0.05 --volt-path a.dat --isi-path isi.txt --ras-path ras.txt

After run the program.

The file `a.dat` will contain voltage traces. In raw `double` binary format. The order is:<a name="V-order"></a>

		V_1(t = 0)     V_2(t = 0)    ... V_n(t = 0)
		V_1(t = dt)    V_2(t = dt)   ... V_n(t = dt)
		...
		V_1(t = L*dt)  V_2(t = L*dt) ... V_n(t = L*dt)

	where L = floor(T/dt), n is number of neurons, T is simulation time, dt  is output time step (--stv).

The file `isi.txt` contains mean Inter-Spike-Intervals for each neurons. In text format.

The file `ras.txt` contains spike events (of the neurons inside the network). In text format. The first column shows neuron indexes (range from 1 to n), second column shows the timings.


### Using the GNU Octave/Matlab interface `gen_neu.m`

(Note: compile the utility function first! See [Build the Matlab interface](#build-matlab))

The interface `gen_neu.m` is a wrapper for `bin/gen_neu`.

```matlab
	pm = [];
	pm.neuron_model = 'HH-PT-GH';  % e.g. LIF-G, LIF-GH, HH-G, HH-GH.
	                       % See bin/gen_neu --help for a complete list.
	pm.net  = [0 1; 0 0];  % Can also be a text file path.
	pm.nI   = 0;           % default: 0. Number of Inhibitory neurons.
			               %             Indexes are later half
	pm.scee_mV = 0.05;
	pm.scie_mV = 0.00;       % default: 0. Strength from Ex. to In.
	pm.scei_mV = 0.00;       % default: 0. Strength from In. to Ex.
	pm.scii_mV = 0.00;       % default: 0.
	pm.pr      = 1.6;        % Poisson input rate.
	pm.ps_mV   = 0.04;       % Poisson input strength.
	pm.t    = 1e4;        % Simulation time.
	pm.dt   = 1.0/32;     % Simulation time step. Default: 1/32.
	pm.stv  = 0.5;        % Output time step, must be multiples of dt.
	pm.seed = 'auto';     % default: 'auto'(or []). Also accept integers.
	pm.extra_cmd = '';    % Optional: put all other command line options here.
	[V, ISI, ras] = gen_neu(pm, 'rm');
```

Now you got voltage traces `V`, mean Inter-Spike-Intervals `ISI` and spike events `ras`.

See [previous section](#V-order) for the format of `V`, `ISI` and `ras`.
Note that there is a transpose for the format of `V` here.

The parameter `'rm'` for `gen_neu` means delete the temporary files after the result been read (just before function return). Otherwise, in next run, with the same parameters, `gen_neu` will use the cached results.


Simualtor Speed
---------------

Time complexity for the simulators:

Simulator     | Number of calls to "Quiet Step"
:-------------|:--------------------------------------
IF-jump       | T * n * (1/dt + pr + n * sp * fr)
simple        | T * n * (1/dt + pr)
SSC           | T * n * (1/dt + pr + (2 + (0~1) * pr*dt) * n * fr)
SSC-Sparse    | T * n * (1/dt + pr + (2 + (0~1) * pr*dt) * n * sp * fr)
SSC-Sparse2   | same as 'SSC-Sparse'
big-delay     | T * n * (1/dt + pr + n * sp * fr)
big-net-delay | same as 'big-delay'

Notes:
* T : simulation time.
* n : number of neurons.
* dt: time step, a fixed value.
* pr: poisson input rate.
* sp: network sparsity.
* fr: neuron mean firing rate.
* "Quiet Step": roughly a RK4 step, it means a call to neuron_model.NextStepSingleNeuronQuiet(), which may contain zero (when the neuron is fully in refractory period), one or two (spiked and awaked from refractory in one time step) calls to RK4 procedure.
* Number of synaptic interaction is: T * n * n * sp * fr. Therefore the time cost can no smaller than this as n goes very large, for all listed simulators.
* For SSC type simulators, if pr * dt * n * sp * fr * dt >> 1, then smaller dt can make the simulation faster. For the SSC-Sparse variation, replace n to n * sp.
* The factor (0~1) in front of pr*dt relates to T, dt, fr and dynamics. For a very rough estimation, it is about 1 - n * fr * dt / 4. The more synchronous the network, the smaller the factor.

### Single neuron time test.

Measured through the interface `gen_neu.m`. See `test/speed_benchmark*` scripts, and [test/speed_benchmark.txt](test/speed_benchmark.txt).
(commit c4562786)

       n = 1+0, t = 100 s, dt = 1/32 ms, stv = 0.5 ms

Model \ sec \ pr      |  1    |  8    |  64   |  512
:---------------------|:-----:|:-----:|:-----:|:-----:
IF-jump               | 0.261 | 0.337 | 0.884 | 5.499
LIF-G + SSC           | 0.619 | 0.714 | 1.424 | 6.672
LIF-GH + SSC          | 0.740 | 0.828 | 1.685 | 8.428
HH-G + SSC            | 1.507 | 1.772 | 3.936 | 21.06
HH-GH + SSC           | 1.574 | 1.857 | 4.140 | 22.09
HH-PT-GH + SSC        | 1.678 | 1.973 | 4.243 | 22.14
HH-GH-sine + SSC      | 1.941 | 2.240 | 4.617 | 23.26
HH-GH-cont-syn        | 8.197 | 9.862 | 23.51 | 132.0

Note: ps_mV = 0.9 / pr.

For comparing `simple` and `SSC` simulator.

Model \ sec \ pr      |  1    |  8    |  64   |  512
:---------------------|:-----:|:-----:|:-----:|:-----:
LIF-G + simple        | 0.360 | 0.440 | 1.115 | 6.361
LIF-GH + simple       | 0.441 | 0.553 | 1.408 | 8.062
HH-G + simple         | 1.181 | 1.453 | 3.604 | 20.56
HH-GH + simple        | 1.258 | 1.547 | 3.816 | 21.89
HH-PT-GH + simple     | 1.487 | 1.705 | 4.117 | 22.76
HH-GH-sine + simple   | 1.644 | 1.925 | 4.285 | 22.85

Legancy simulators ([raster_tuning](https://bitbucket.org/bewantbe/ifsimu) based).

Model \ sec \ pr       |  1    |  8    |  64   |  512
:----------------------|:-----:|:-----:|:-----:|:-----:
legancy-LIF-G          | 0.797 | 0.920 | 1.790 | 7.355
legancy-LIF-GH         | 1.000 | 1.157 | 2.326 | 10.13
legancy-HH-GH-cont-syn | 2.254 | 2.603 | 5.413 | 26.63

### Network

  1000 neurons, fully connected.

  n = 800+200, t = 0.1 s, dt = 1/32 ms, stv = 0.5 ms

model                  | sec   | mean freq (Hz)
:----------------------|:-----:|:------:
IF-jump                | 0.458 | 30.290
LIF-G + SSC            | 0.818 | 24.430
LIF-GH + SSC           | 1.084 | 24.480
HH-G + SSC             | 3.692 | 34.570
HH-GH + SSC            | 4.041 | 34.720
HH-PT-GH + SSC         | 4.068 | 34.490
HH-GH-sine + SSC       | 5.354 | 34.720
||
LIF-GH + simple        | 0.679 | 24.480
HH-GH + simple         | 1.599 | 34.700
||
legancy-LIF-G          | 1.541 | 39.660
legancy-LIF-GH         | 1.786 | 34.690
legancy-HH-GH-cont-syn | 244.0 | 37.380


  1000 neurons, 10% connected.

  n=800+200, t = 0.1 s, dt = 1/32 ms, stv = 0.5 ms

model                  | sec   | mean freq (Hz)
:----------------------|:-----:|:------:
IF-jump                | 0.339 | 25.820
||
LIF-GH + simple        | 0.571 | 26.580
LIF-GH + SSC-Sparse    | 0.657 | 26.580
LIF-GH + SSC-Sparse2   | 0.659 | 26.580
LIF-GH + SSC           | 1.052 | 26.580
||
HH-GH + simple         | 1.473 | 35.970
HH-GH + SSC-Sparse     | 1.808 | 35.970
HH-GH + SSC-Sparse2    | 1.749 | 35.970
HH-GH + SSC            | 4.019 | 35.970

  Network of different size. Sparsity is 10%.
  Firing rates are around ~30 Hz, like above.

model \ sec \ n       |  100   |  300  | 1000  | 3000  | 10000
:---------------------|:------:|:-----:|:-----:|:-----:|:-----:
IF-jump               |  0.029 | 0.063 | 0.293 | 1.840 | 17.865
LIF-GH + simple       |  0.092 | 0.156 | 0.515 | 2.349 | 18.543
LIF-GH + SSC-Sparse   |  0.095 | 0.168 | 0.592 | 2.940 | 24.801
HH-GH + simple        |  0.198 | 0.427 | 1.391 | 4.947 | 27.095
HH-GH + SSC-Sparse    |  0.198 | 0.462 | 1.674 | 7.415 | 55.261
LIF-GH + SSC          |  0.101 | 0.210 | 1.043 | 6.748 | 69.236
HH-GH + SSC           |  0.226 | 0.665 | 4.047 | 27.71 | 284.15

