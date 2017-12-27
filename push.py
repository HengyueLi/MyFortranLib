#!/usr/bin/env python3

import datetime,os



time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
gitadd = "git add .;"
gitcom = "git commit -m "+"\""+time+"\";"
gitpus = "git push -u origin master --force;"

os.system(  gitadd + gitcom + gitpus )
