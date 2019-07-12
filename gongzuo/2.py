#!/uar/bin/python3
import subprocess
def runcmd(command):
    ret=subprocess.run(command,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,encoding="UTF-8",timeout=1)
    if ret.returncode==0:
        print("success:",ret)
    else:
        print("error:",ret)
runcmd(["dir","/b"])
runcmd(" exit 1")