cd/
mkdir datapond
echo alias "lt='ls -lhtr'" > ~/.bashrc && source ~/.bashrc
cd ~/datapond/Dependencies
wget -nc https://repo.anaconda.com/miniconda/Miniconda3-py311_23.5.0-3-Linux-x86_64.sh 
sudo bash Miniconda3-py311_23.5.0-3-Linux-x86_64.sh -b -p ~/datapond/Dependencies/conda
echo 'export PATH=$PATH:/home/ubuntu/datapond/Dependencies/conda/bin' >> ~/.bashrc
source ~/.bashrc
sudo chown -R $USER /home/ubuntu/datapond/Dependencies/conda

conda install -y python=3.8.12
sudo apt-get install -y cutadapt cmake htop libharfbuzz-dev libfribidi-dev libcurl4-openssl-dev libssl-dev libfontconfig1-dev libproj-dev s3fs
cd ~/Scripts
conda env create -n bridge_env -f bridge_conda_env.yml
# create reqirements.txt from the environment
conda init
source ~/.bashrc 
conda activate bridge_env
curl -s https://get.nextflow.io | bash
echo 'export PATH=$PATH:/home/ubuntu/Scripts' >> ~/.bashrc
source ~/.bashrc 

wget https://github.com/alexdobin/STAR/archive/2.7.11b.tar.gz
tar -xzf 2.7.11b.tar.gz
cd STAR-2.7.11b

cd STAR/source
make STAR
echo 'export PATH=$PATH:/home/ubuntu/Scripts/STAR-2.7.11b/source' >> ~/.bashrc


