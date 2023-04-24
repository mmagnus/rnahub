Step 1: Install software programs

1) Install nhmmer, which is included in the HMMER software package. Download the tar archive from http://hmmer.org/
Per the User Manual, install using the following commands. 

	tar xf hmmer.tar.gz
	cd hmmer-3.3.2
	./configure
	make
	make install

You will also need to install the Easel commands as part of your installation, using this:

	cd easel; make install

2) Install R-scape. Download the tar archive from http://rivaslab.org/software.html . Per the User Manual, 
install using the following commands:

	tar xf rscape.tar.gz
	cd rscape_v1.5.16
	./configure
	make
	make install
