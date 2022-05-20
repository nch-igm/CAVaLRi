from urllib.parse import urlparse

import os, shutil, boto3
from awscli.clidriver import create_clidriver


def validate_s3_path(s3_url):
    """Validate an s3 path."""
    if not s3_url.startswith("s3://"):
        raise ValueError("Invalid s3 path.")
    parse = urlparse(s3_url)
    bucket, path = parse.netloc, parse.path.lstrip('/')
    return bucket, path


def download_s3_file(s3_url, local_file):
    """Download file from s3 locally."""
    bucket, key = validate_s3_path(s3_url)
    client = boto3.client("s3")
    client.download_file(bucket, key, local_file)


def upload_s3_file(local_file, s3_url):
    """Upload a local file to s3."""
    bucket, key = validate_s3_path(s3_url)
    client = boto3.client("s3")
    client.upload_file(local_file, bucket, key)


def download_s3_directory(s3_dir_name, local_dir_name):
    """Validate the s3 url, create cmd, download."""
    validate_s3_path(s3_dir_name)
    cmd = ['s3', 'cp', '--recursive', s3_dir_name, local_dir_name]
    rc = create_clidriver().main(cmd)
    if rc != 0:
        raise RuntimeError(f"There was an error downloading the directory. {rc}")


def upload_s3_directory(local_dir_name, s3_dir_name):
    """Validate the s3 url, create cmd, upload."""
    validate_s3_path(s3_dir_name)
    cmd = ['s3', 'cp', '--recursive', local_dir_name, s3_dir_name]
    rc = create_clidriver().main(cmd)
    if rc != 0:
        raise RuntimeError(f"There was an error uploading the directory. {rc}")

 
def flex_input(in_path: str, out_dir: str=".", directory: bool=False):
    """
    Allows specifying a path to either an S3 file or a local file.
    If valid local file, returns the same path as input.
    If S3 file, downloads the file and returns the new local path.
    If input is a directory, must set 'directory' parameter to True.
    """
    if in_path is None:
        out_path = in_path
    else:
        if directory is True:
            if in_path.startswith("s3://"):
                out_path = f'{out_dir}/{os.path.basename(in_path)}'
                print(out_path)
                download_s3_directory(in_path,out_path)
            else:
                if os.path.isdir(in_path):
                    out_path = in_path
                else:
                    raise ValueError(f"Path {in_path} does not specify either a valid S3 path or a valid local directory.")
        else:
            if in_path.startswith("s3://"):
                out_path = f'{out_dir}/{os.path.basename(in_path)}'
                print(out_path)
                download_s3_file(in_path,out_path)
            else:
                if os.path.isfile(in_path):
                    out_path = in_path
                else:
                    raise ValueError(f"Path {in_path} does not specify either a valid S3 path or a valid local file.")
    return out_path
    

def flex_output(in_path: str, out_dir: str=".", file_rename: str = None):
    """
    Copies a file or directory to either a local output directory or an S3 path.
    If valid local path, returns the same path as input.
    If S3 path, uploads the input file or directory to that path.
    """
    
    if file_rename is None:
        out_path = f'{out_dir}/{os.path.basename(in_path)}'
    else:
        out_path = f'{out_dir}/{file_rename}'
    
    if in_path == out_path:
        print("Source and destination are the same. No copy needed.")
        
    elif out_dir.startswith("s3://"):
        print(out_path)
        if os.path.isdir(in_path):
            upload_s3_directory(in_path, out_path)
        elif os.path.isfile(in_path):
            upload_s3_file(in_path, out_path)
        else:
            raise ValueError(f"Input {in_path} is not a valid file or directory path.")
    else:
        if os.path.isdir(out_dir):
            if os.path.isdir(in_path):
                shutil.copytree(in_path,out_path)
            else:
                shutil.copy(in_path, out_path)
        else:
            raise ValueError(f"Output {out_path} does not specify either a valid S3 path or a valid local file.")
    return out_path