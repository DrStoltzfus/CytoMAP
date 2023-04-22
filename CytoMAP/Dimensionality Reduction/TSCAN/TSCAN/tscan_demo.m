data = csvread('dataframe.csv',1,1);
procdata = preprocess(data);
lpsmclust = exprmclust(procdata);    
plotmclust(lpsmclust);


    