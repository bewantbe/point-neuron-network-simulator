Accurate Point Neuronal Network Simulator
=========================================
(Note: **NOT** for artifical neural network)

This program aims to provide an accurate simulator and a portable code framework for point neuron model networks.

A handy GNU Octave/Matlab interface to the core simulator is also provided.

The accuracy is achieved by breaking the simulation time steps (dt) according to the (e.g. spike) event timings. So that the dynamics inside the sub-time-steps are smooth and the classical ODE solves (notably RK4) work well.

Currently there are (parenthesis: Name used inside the code):

  * Leaky Integrate-and-Fire models with exponential decay conductance (LIF-G).
  * Leaky Integrate-and-Fire models with "alpha function" conductance (LIF-GH).
  * Hodgkin Huxley neuron model with exponential decay conductance (HH-G).
  * Hodgkin Huxley neuron model with "alpha function" conductance (HH-GH).

Optional:

  * Add an alternating current to these models.
  * Add synaptic delay. (Currently constant delay, and must larger than dt)

See `doc/neuron_models.pdf` for model details.

You can provide external input events in a file or just let the program generate (Poisson) for you.

Comparison to other simulators
------------------------------

* NEST

    NEST is a full featured point neuron simulator. While this one (APNNS) has much less neuron models and tunable settings. The function of APNNS is extended by directly modifying the source code.

    NEST is a big project. APNNS is rather small, roughly about 3 thousand lines of C++ code.
    
    NEST use a high order ODE solver for each neurons (by default it is Runge–Kutta–Fehlberg method (aka rkf45 or ode45()) through GNU Scientific Library). But the interactions between neurons has only order one accuracy. Because the interactions is performed at the boundary of time step. While APNNS will deal with this case accurately.

* NEURON

    NEURON is a more detailed/physiology simulator. APNNS only simulate "point" neurons.

* Brian

    (Not yet compared)


Build from source
-----------------

### Build under linux

You will need to install libraries Eigen and Boost.ProgramOptions first.
Then simply `make`, you will get an executable `bin/gen_neu`.

There are matlab scripts in `mfile/`, notably the interface `mfile/gen_neu.m`.
To use the matlab interface, you need to compile `mfile/BKDRHash.c` and `mfile/randMT19937.cpp`. See the source code for how to compile them.


### 在 Windows 下用 MSVC 编译

示例环境：

  Win 10, Microsoft Visual Studio Community 2015

步骤：

1. 下载 boost

		http://www.boost.org/users/download/
		点 Prebuilt windows binaries.
		(https://sourceforge.net/projects/boost/files/boost-binaries/1.60.0/)
		VS -> Help -> About 查看 VS 的版本。 比如 2015 是 version 14.
		下载(32bit 通用性好一些) boost_1_60_0-msvc-14.0-32.exe
		安装，按默认路径。

2. 下载 Eigen

		http://eigen.tuxfamily.org/
		选择 "latest stable release" 
		（http://bitbucket.org/eigen/eigen/get/3.2.8.zip）
    
		解压复制到编译系统文件夹
			C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\include
		重命名文件夹为 eigen3
		(也就是应该有文件 C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\include\eigen3\Eigen\Dense)
	
		也可以放到任意放一个地方（路径最好不要含中文）, 然后添加 Include Directories

3. 用 MSVC 打开 point-neuron-network-simulator.sln
  或手动新建工程，添加需要的文件。
  
  * 设置 boost

			如非默认路径，需手动在工程中添加 include 目录：
			在工程属性里设定
				Include Directories 为 C:\local\boost_1_60_0
				Library Directories 为 C:\local\boost_1_60_0\lib32-msvc-14.0
			如果是 64bit 工程，则
				Library Directories 为 C:\local\boost_1_60_0\lib64-msvc-14.0

  * 设置 Eigen

		如库文件没有放到编译系统文件夹，请自行设定 Include Directories。

	测试用命令行参数

		--neuron-model HH-GH --net - --dt 0.03125 --t 1e3 --pr 1 --ps 0.05  --seed 1

	没报错就是编译正常了。

