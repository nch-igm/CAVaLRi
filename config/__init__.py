import os
import yaml

# def get_conf(config_path_list=[
#                     os.path.join('..', '..', 'config', 'config.yaml'),
#                     os.path.join('..', 'config', 'config.yaml'), 
#                     os.path.join('config', 'config.yaml')
#                     ]):

def get_conf():

    # Initiate configuration object
    config = {}

    root_path = os.path.abspath(os.path.dirname(__file__))
    config_path_list = [os.path.join(root_path, 'config.yaml')]

    for cf in config_path_list:
        try:
            with open(cf,'r') as f:
                conf = yaml.load(f, Loader=yaml.FullLoader)
                config.update(conf)
        except:
            pass

    return config

config = get_conf()
