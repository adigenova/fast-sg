# Helper script to install ScaffMath and their dependencies


echo "ScaffMatch will be installed in $PWD"
INSTALL_DIR=$PWD
BOWTIE2=""
PYTHON2=""
PYTHON2B=""
PYTHONPATH2=""

echo "Downloading ScaffMatch"
wget -c http://alan.cs.gsu.edu/scaffmatch/ScaffMatch-0.9.tar.gz
echo "Checking bowtie 2"
if command -v bowtie2 >/dev/null 2>&1 ; then
	echo "bowtie2 found"
	echo "version: $(bowtie2 --version)"
	BOWTIE2=$(dirname $(command -v bowtie2))
else
 	wget -c https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.4/bowtie2-2.3.4-linux-x86_64.zip/download -O bowtie2-2.3.4-linux-x86_64.zip
	unzip bowtie2-2.3.4-linux-x86_64.zip
	BOWTIE2=$INSTALL_DIR/bowtie2-2.3.4-linux-x86_64
fi

echo "Checking Python 2.7"
if command -v python2.7 >/dev/null 2>&1 ; then
	echo "python2.7 found"
	echo "version: $(python2.7 --version)"
	PYTHON2=$(dirname $(command -v python2.7))
	PYTHON2B= $(command -v python2.7)
else
	wget -c https://www.python.org/ftp/python/2.7.13/Python-2.7.13.tgz -O Python-2.7.13.tgz
	tar zxvf Python-2.7.13.tgz
	cd Python-2.7.13
	./configure --prefix=$PWD/LOCALPYTHON
	make all
	make install
	PYTHON2=$PWD/LOCALPYTHON/bin/
	PYTHON2B=$PWD/LOCALPYTHON/bin/python2.7
fi

echo "Intalling pip for Python dependencies locally"
cd $INSTALL_DIR
wget https://bootstrap.pypa.io/get-pip.py && $INSTALL_DIR/Python-2.7.13/LOCALPYTHON/bin/python get-pip.py --user

echo -n "Checking for Networkx python library: "

$PYTHON2B -c "import networkx" &> /dev/null
if [ $? -eq 1 ]; then
	echo "Installing networkx"
	$HOME/.local/bin/pip install networkx --user
	PYTHONPATH=~/.local/lib/python2.7/site-packages
else
    echo "OK"
fi

echo -n "Checking for Numpy python library: "
$PYTHON2B -c "import numpy" &> /dev/null
if [ $? -eq 1 ]; then
	echo "Installing numpy"
	$HOME/.local/bin/pip install numpy --user
	PYTHONPATH2=~/.local/lib/python2.7/site-packages
else
    echo "OK"
fi


echo -n  "SET PATH to:"
echo -n "export PATH=${BOWTIE2}:${PYTHON2}:$PATH"
echo -n "SET PYTHONPATH to "
echo -n "export PYTHONPATH=${PYTHONPATH2}:${PYTHONPATH}"

