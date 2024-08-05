# Folders in the prac
There has been a bit of confusion of where is what so I thought I put some info on here for you.

## Explanation
There are two different areas on the servers that I saw students work in and this caused issues with the commands that I posted finding things. If you cannot find a file you need to:

1. check where you copied your files
2. check where your command is looking

If you decided to arrange your files differently to the prac instructions (e.g. you added an additional directory to sort it all in) then you have to adjust some of the commands. *chatGPT is fab at picking these things*

You need to distinguish between files that we provide (don't change those paths) and files that you have uploaded or created (these are the ones that the command then complains about not finding).

The last thing was that we sometimes set a directory variable `sourcedir`. These get lost when you log out and in or if you open a new terminal. So be sure to rerun this `sourcedir = /path/to/whatever/you/need`

## Areas on server
### The log on area - home
When you log onto the server with ssh you tend to see the home directory which is your home directory. See it as your landing folder that can pollute quite quickly as everything defaults to wanting to send things there.

The path looks like this if you have account s-99: `/home/s-99/`

And when you use the `cd ~` it will take you there.

**Don't put your files here but check if a file went here if it copied and you cannot see it in your account!**

### The account folder
You may have noticed that your files and other things were always copied and created into your *account* folder. This is were you can put all the files that you have accumulated. So copy your files here (e.g. `Av.chloe.gff` etc)

This folder should also contain the `chloe` folder that then the julia command piles in those gff3 files.

It looks like this if you have account s-99: `/mnt/s-ws/s-99`

(mnt stands for mount so this is a data volume that is better to "pollute" with files.)


### Shared area
Is also on the mount and this is where our source directory points to or where we have some files that you need to copy across or access.

An example is `/mnt/s-ws/everyone`


## File transfer and File finding
I think the easiest is to you VS code to transfer files. You can do all tasks in VS code:

### Connect to server
In VScode find the green arrows button on the bottom left and choose `connect to host` (in current window if you don't want to open a new one):

![connect to server](pic/vs1.png)

Then enter your account and your server. You will get prompted for you password.

![connect to server](pic/vs2.png)

### Open your directory to not get lost in home
To get to your directory you click on the blue button on the left `open folder`. It will open the dialog and you will see the default is e.g. `/home/s-99`or whatever you account is.

![connect to folder](pic/vs3.png)

This is where you can set your folder to the account area by giving it the right path. Then click ok

![connect to folder](pic/sv4.png)

You can see now on the left the folder contents. Here you can right click on a file and download it without using the command. In order to see another area you have to open a new folder view.

### Filezilla
If you are copying files no matter what you use you need to make sure where your app is pointing at. In Filezilla or alike apps you need to check where it is putting your files. Check the path given on the top and if it sais `/home/s-99` you will transfer it to home not the account folder.

**Make sure your files are all in the account area** 
