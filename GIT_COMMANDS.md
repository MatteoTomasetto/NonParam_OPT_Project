## Clone the github directory
cd "C:/Users/matte/Desktop"  # set the path where we want to add the directory

git clone https://github.com/teotom10/NonParam_OPT_Project.git

## Commands to update the files

cd "C:/Users/matte/Desktop/NonParam_OPT_Project"  # select the wanted directory

git pull  # update the folder and files on own PC

git add -A                       # add all files on remote folder
git add "file_name.extension"    # add one file on remote folder
git commit -m "Commit_text"      # commit a change
git push origin head:master      # update the remote folder on github 

## ITER TO FOLLOW
# 1) Modify and work on local file(s)
# 2) git pull
# 3) Update the file(s) in the shared folder
# 4) git add + git commit + git push