Major updates on January 2, 2013

1) Ability to read metamorph file as well as loading ND file
2) Revised commandline inputs (see below) 
3) Added MaskGenerating module for curating mask of each cell

//How to use modules//

CellTracking (or FRETPlotting/MaskGenerating)
varargin{1} - filetype 1 = PE export, 2 = tifstack, 3 = metamorph images
varargin{2} - trackinginfo  - [templateCH tp_1 tp_end]       
varargin{3} - imagelocation - [row col field plane]
varargin{4} - channelnames eg. {'CFP';'FRET';'RFP'};
varargin{5} - fileformat or tiffstack file name  
examples: 
CellTracking(1,[templateCH tp_1 tp_end],[row col field plane])
CellTracking(2,[templateCH tp_1 tp_end],[],{'1';'2';'3'},<tiff stack file>)
CellTracking(3,[templateCH tp_1 tp_end],[],{'CFP';'YFP';'RFP'},'2012-11-16_%s_xy087_t%03g.tif')

Special command for Metamorph input:

CellTracking(3,'12322012.nd')


//Archives of changes//

Major updates on December 2, 2012

1) Ability to read image file (file type 3)
2) Add 'Editable' button. Track points can only be edited when this toggle switch is pressed.
3) Faster tracking and playing mode (cutting down process time from not initiating track points).
4) Revised commandline inputs (see below) 
5) Automatic detection of image display range. Threshold can be changed by editing min and max threshold.