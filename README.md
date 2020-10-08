# What's MetaCache-MPI about? #

**MetaCache-MPI** is a software tool that allows to use the metagenomics minhashing algorithm from [**metacache**][1] by means of a High Performance Computing environment by using the MPI library.

If you use **MetaCache-MPI**, please cite this article:

> José M. Abuín, Nuno Lopes, Luís Ferreira, Tomás F. Pena and Bertil Schmidt. ["Big Data in metagenomics: Apache Spark vs MPI"][4]. PLoS ONE 15(10). 06 October 2020.

The original **MetaCacheSpark** and **metacache** articles are:

> Robin Kobus, José M. Abuín, André Müller, Sören Lukas Hellmann, Juan C. Pichel, Tomás F. Pena, Andreas Hildebrandt, Thomas Hankeln and Bertil Schmidt. ["A Big Data Approach to Metagenomics for All-Food-Sequencing"][3]. BMC Bioinformatics, 21, Article number: 102 (2020).

> André Müller, Christian Hundt, Andreas Hildebrandt, Thomas Hankeln, Bertil Schmidt. ["MetaCache: context-aware classification of metagenomic reads using minhashing"][2]. Bioinformatics, Volume 33, Issue 23, 01 December 2017, Pages 3740–3748.

# Getting started #

## Requirements
In order to build and run **MetaCache-MPI** the following items are needed:

* CMake >= 3.10.
* MPI. (Tested with OpenMPI 3.0.2 and 2.1.6).
* The Pthreads library.
* A C++ compiler able to build C++14 code. (Tested with G++ 7.3.1 and 6.4.0)

## Building
The default way to build **MetaCache-MPI** is:

	git clone https://github.com/jmabuin/metacache-mpi.git
	cd metacache-mpi
	mkdir build && cd build
   	cmake -DMEDIUM_TARGETS=True ..
   	make

This will create the *metacache_mpi* file.

## Launching
Depending on the queue system installed in the cluster we are using, the launching commands can be slighlty different. Examples of how to launch **metacache-mpi** for *building* and *querying* (the only two implemented commands) can be found in the *script* folder.

Regarding the parameters related to the algorithm itself, **metacache-mpi** uses are the same than **metacache**. Information about these parameters can be found at [https://github.com/muellan/metacache][1]


[1]: https://github.com/muellan/metacache
[2]: https://doi.org/10.1093/bioinformatics/btx520
[3]: http://dx.doi.org/10.1186%2Fs12859-020-3429-6
[4]: https://doi.org/10.1371/journal.pone.0239741