import subprocess

result = subprocess.run(["/home/ubuntu/DE_NOVO_PROTEINS/SCRIPTS/martini_shell.sh"], capture_output=True, text=True)

print("stdout", result.stdout)
print("stderr", result.stderr)

#process = subprocess.Popen(["/home/ubuntu/DE_NOVO_PROTEINS/SCRIPTS/martini_shell.sh"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#output, error = process.communicate()
#print('output', output)
#print('error', error)

