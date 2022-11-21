import os
import shlex
import subprocess
import argparse
from config import *

def worker(cmd, **kwargs):
    if 'shell' in kwargs.keys():
        if kwargs['shell']:
            p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell = True)
    else:
        cmd = shlex.split(cmd)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p.wait()
    out, err = p.communicate()
    print(out.decode() if out else err.decode())
    return out.decode() if out else err.decode()

def main(dir):
    
    # Set to default, transfer to absoute path
    dir = os.path.abspath('dependencies') if not dir else os.path.abspath(dir)
    
    # Create the directory if it doesn't exist
    if not os.path.exists(dir):
        os.mkdir(dir)
    
    # Download bcftools
    print('Installing bcftools')
    link = 'https://github.com/samtools/bcftools/releases/download/1.16/bcftools-1.16.tar.bz2'
    bcftools_dir = os.path.basename(link)
    bcftools_dir = os.path.join(dir, bcftools_dir[:bcftools_dir.find('.tar')])
    worker(f'wget {link} -P {dir}')
    worker(f"tar -xjf {os.path.join(dir, 'bcftools-1.16.tar.bz2')} -C {dir}")
    worker(f"""
        cd {bcftools_dir};
        ./configure --prefix={dir};
        make;
        make install
    """, shell = True)

    # Download htslib
    print('Installing htslib')
    link = 'https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2'
    htslib_dir = os.path.basename(link)
    htslib_dir = os.path.join(dir, htslib_dir[:htslib_dir.find('.tar')])
    worker(f'wget {link} -P {dir}')
    worker(f"tar -xjf {os.path.join(dir, 'htslib-1.16.tar.bz2')} -C {dir}")
    worker(f"""
        cd {htslib_dir};
        ./configure --prefix={dir};
        make;
        make install
    """, shell = True)
    
    # Download annovar databases
    if config['genome_build'] == 'hg38':
        dbs = ['refGene','gnomad30_genome','clinvar_20220320']
        for db in dbs:
            print(f'Downloading annovar db: {db}')
            cmd = f"{os.path.join(config['annovar_scripts'], 'annotate_variation.pl')} -downdb -buildver hg38 -webfrom annovar {db} {os.path.join(dir, 'human_db')}"
            # worker(cmd)

    # Download LIRICAL
    # link = 'https://github.com/TheJacksonLaboratory/LIRICAL/releases/download/v1.3.4/LIRICAL.jar'
    # print('Downloading LIRICAL.jar')
    # worker(f'wget {link} -O -P {dir}')
    print('Downloading LIRICAL data')
    worker(f"java -jar {os.path.join(dir, 'LIRICAL.jar')} download -d {os.path.join(dir, 'lirical_data')}")
    
    # Download Exomiser
    print('Downloading Exomiser data')
    link = 'https://data.monarchinitiative.org/exomiser/latest/2209_hg38.zip'
    worker(f'wget {link} -P {dir}')
    

if __name__ == "__main__":
    
    # Intialize parser
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--dir', '-d', type=str, help='Path where dependencies will be downloaded to. If not provided, a "depencies" folder will be created in the root directory of the project')
    args = parser.parse_args()

    x = main(args.dir)
