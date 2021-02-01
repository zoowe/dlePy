# Step 1:
Prepare your input. For each calculation, basic input files must present in separate folder.

# Step 2:
Prepare job script file with ```#SBATCH -J $name``` 

# Step 3:
Distribute job script file to each folder. This task can be done with ```cp.py```

# Step 4:
Edit ```job_monitor.py``` as your wish. If you want to have email function working, you need to set the following enviroment variables:
You can use mailjet.com for email service
```
export mail_server='send_email_server'
export mail_port='587'
export mail_login='send_email_username'
export mail_password='send_email_password'
export from_address='youremail'
export to_address='youremail'
```

Note that you can define yourown function to determine when your job is done to replace ```if_vasp_done```.

# Step 5:
Run ```nohup python3 -u job_monitor.py >& job_monitor.log &```
View job_monitor.log for details.

