# TCPDat
This is the repository for the Tropical Cyclone Precipitation Dataset, or TCPDat. 

 TCPDat: Tropical Cyclone Precipitation Dataset
 Public version 1.0 release 1.0
 Date updated: 03/02/2020 (dd/mm/yyyy)

 If you use this code or the pre-generated TCPDat data, please include the
 following citation (or whichever citation style your journal of choice requies):

 Bregy, J.C., J.T. Maxwell, S.M. Robeson, J.T. Ortegren, P.T. Soulé, and
 P.A. Knapp, 2019: Spatiotemporal variability of tropical cyclone
 precipitation using a high-resolution, gridded (0.25° x 0.25°) dataset
 for the eastern United States, 1948–2015. Journal of Climate, 
 33(5), 1803–1819, doi:10.1175/JCLI-D-18-0885.1.

 Stable link to the original paper: https://doi.org/10.1175/JCLI-D-18-0885.1


 This code generates a gridded tropical cyclone precipitation (TCP)
 product for the eastern United States (i.e., east of the Rocky Mountain
 continental divide) from 1948 to present-1 (or present-2, depending on
 the DOY). This script requires two (2) publically available datasets in
 order to create TCPDat. The first is HURDAT2, specifically the North
 Atlantic version of IBTrACS (International Best Track Archive for Climate
 Stewardship). The script is designed to work with the netCDF version of
 IBTrACS. You can find the latest version of IBTrACS at:
 Cite IBTrACS accordingly.

 The other dataset required is the Climate Prediction Center Unified
 Guage-Based Analysis of Daily Precipitation over CONUS (CPC US Unified; URD)
 from the Physical Sciences Division of the Earth Systems Research
 Laboratory. The latest version can be found here: 
 https://www.esrl.noaa.gov/psd/thredds/catalog/Datasets/cpc_us_precip/catalog.html
 (Note: the years 2007 onward are located in the RT subfile. CPC started doing something different.
 It's inconvenient, but such is life.) Like IBTrACS, CPC US Unified should be a netCDF. You will have to
 download gridded precipitation for each year that you are interested in.
 The script is written to cover every year starting at 1948 (per CPC US
 Unified). In future iterations, I plan to have query statements that
 allow you to define your time range, but plans always look great on paper.
 Anyway, cite CPC US Unified accordingly.

 This is my attempt at annotating this code. It is going to be rough; this
 is literally the first time I have ever done something like this. Let the
 GitHub gods have mercy upon my digital soul. If there are any questions
 about anything that is in the code, or if you have any suggestions to
 improve the code, please do not hesitate to get in touch with me. My most
 frequently checked email is jbregy@indiana.edu, but I also have a stable
 address at joshua.bregy@gmail.com. I suppose you could get in touch via
 GitHub, ResearchGate, or Twitter (@prehistormic) if you would like.

 READ THIS BEFORE RUNNING THE SCRIPT
 You will need to create several folders for this script to work. Each
 folder houses different data, and the code moves between different folders 
 to either use or save datasets. I wrote it this way to force myself to be
 organized (normally, I'm not, but I knew my future self would genuinely
 appreciate it). Here are the folders that you will need:

 1) A folder containing your CPC US Unified files. I called my
 CONUS_daily_precip_1948_2016. (I accidentally downloaded p2016.nc, but
 HURDAT2 had not been updated to include 2016.) I will be changing the
 file name in the command to CONUS_daily_precip, so I suggest that you
 create a file with the same name. That, or change the name in the code
 yourself. Alternatively, you could just comment out directory changes
 that are used to access data as long as the folder containing your
 precipitation (and I suppose IBTrACS, too) is included whenever you set
 your path.

 2) A folder for the TCP data. I have a line in the code that creates that
 folder. It is named TCP_data. 
