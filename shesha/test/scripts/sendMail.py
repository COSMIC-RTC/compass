import smtplib
from email.mime.text import MIMEText
smtp = smtplib.SMTP('smtp.obspm.fr')
msg = MIMEText('Simulation finished.')
msg['From'] = 'micmac'
msg['To'] = 'Script runner'
msg['Subject'] =  "Script finished!!!!!"
smtp.sendmail('micmac@obspm.fr', ["fabrice.vidal@obspm.fr"], msg.as_string())
