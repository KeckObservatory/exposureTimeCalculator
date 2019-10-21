#
# Script to manager (start/stop/restart) flask servers
#
import argparse
import os
import subprocess
import psutil
import getpass

def is_server_running(server):
    '''
    Returns PID if server is currently running, else 0
    '''

    current_user = getpass.getuser()
    list1 = ['python', server]
    for proc in psutil.process_iter():
        pinfo = proc.as_dict(attrs=['name', 'username', 'pid', 'cmdline'])
        if pinfo['username'] != current_user:
            continue

        list2 = pinfo['cmdline']
        res = [name for name in list1 if any(name in s for s in list2)]
        if res == list1: return pinfo['pid']

    return 0

def process_stop(pid):
    '''
    Use psutil to kill the process ID
    '''

    if pid == 0:
        print(server, 'is not running')
    else:
        print('Killing PID', pid)
        p = psutil.Process(pid)
        p.terminate()
        pid = 0

    return pid

def process_start(pid, server):
    '''
    Start the requested server
    '''

    if pid > 0:
        print(server, 'already running with PID', pid)
    else:
        print('Starting', server)
        cmd = ['/usr/local/anaconda3-5.0.0.1/bin/python', server]
        p = subprocess.Popen(cmd)


# Define input parameters

parser = argparse.ArgumentParser(description='manager.py input parameters')
parser.add_argument('server', type=str, help='flask server name')
parser.add_argument('command', type=str, help='start, stop or restart')

# Get input parameters

args = parser.parse_args()
server = args.server
command = args.command

# Verify command

assert command in ['start', 'stop', 'restart'], print('Incorrect command')

# Determine directory

dir = os.path.dirname(os.path.realpath(__file__))
os.chdir(dir)

# Check if server file exists

server = ''.join((server, '.py'))
assert os.path.isfile(server), print('server does not exist')

# Check if server is running

pid = is_server_running(server)

# Do the request

if command == 'stop':
    pid = process_stop(pid)
elif command == 'start':
    process_start(pid, server)
elif command == 'restart':
    pid = process_stop(pid)
    process_start(pid, server)

exit()
