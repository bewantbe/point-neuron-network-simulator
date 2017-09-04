Accurate Point Neuronal Network Simulator
=========================================
(Note: **NOT** for artificial neural network)

This program aims to provide an accurate simulator and a portable code framework for point spiking neuron model networks.

A handy GNU Octave/Matlab interface to the core simulator is also provided.

The accuracy is achieved by breaking the simulation time steps (dt) according to the (e.g. spike) event timings. So that the dynamics inside the sub-time-steps are smooth and the classical ODE solves (notably RK4) work well.

Currently there are (parenthesis: Name used inside the code):

  * Leaky Integrate-and-Fire models with exponential decay conductance (LIF-G).
  * Leaky Integrate-and-Fire models with "alpha function" conductance (LIF-GH).
  * Hodgkin Huxley neuron model with exponential decay conductance (HH-G).
  * Hodgkin Huxley neuron model with "alpha function" conductance (HH-GH).
  * Hodgkin Huxley neuron model with continous (but mostly around spike) synaptic interaction (HH-GH-cont-syn).

Optional:

  * Add an alternating current to these models.
  * Add an constant current to these models.
  * Add synaptic delay. (The delay must larger than dt)
  * For HH type neuron, option to use the peak as spike time.

See `doc/neuron_models.pdf` for model details.

You can provide external input events in a file or just let the program generate (Poisson input) for you. Input events can be excititory or inhibitory, and the strength can be specified at each timing.

Refer to `bin/gen_neu --help` for the command line option help (after compilation).

### Note:

   Code still in developing stage, use with care.


Comparison to other simulators
------------------------------

* [NEST](http://www.nest-simulator.org/)

    NEST is a full featured point neuron simulator. While this one (APNNS) has much less neuron models and tunable settings. The function of APNNS is extended by directly modifying the source code.

    NEST is a big project. APNNS is rather small, roughly about 3 thousand lines of C++ code.
    
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
		重命名文件夹为 eigen3
		(也就是应该有文件 C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\include\eigen3\Eigen\Dense)
	
		也可以放到任意放一个地方（路径最好不要含中文）, 然后添加 Include Directories

3. 用 MSVC 打开 point-neuron-network-simulator.sln
  或手动新建工程，添加需要的文件。
  
  * 设置 boost

			如非默认路径，需手动在工程中添加 include 目录：
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

