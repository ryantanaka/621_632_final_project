import subprocess

implementation = "BP"
run_script = "./run.sh {} {}"

print("RUNNING IMPLEMENTATION: {}".format(implementation))

for i in range(4, 24 + 1):
    subprocess.Popen(run_script.format(implementation, 2**i), shell=True).wait()
