library(magick)
library(dplyr)
library(bmp)
library(pixmap)


#read the image
pika<-image_read("http://pngimg.com/uploads/pokemon/pokemon_PNG9.png")
#charmer<-image_read("https://pngimg.com/uploads/pokemon/pokemon_PNG154.png")
#image_info(pika)
#print(pika)
#print(charmer)
#resize the image
#recommended size is 100 ~ 200
pika<-image_resize(pika,"150x")
#charmer<-image_resize(charmer,"150x")
#convert to bmp format
pika_bmp<-image_convert(pika,"bmp")
#image_info(pika_bmp)
image_write(pika_bmp,path="/Users/cherry/Desktop/pika.bmp",format="bmp")

#charmer_bmp<-image_convert(charmer,"bmp")
#image_info(charmer_bmp)
#image_write(charmer_bmp,path="/Users/cherry/Desktop/charmer.bmp",format="bmp")

#convert to black and white image
#by PhotoScape X

dot.image<-function(file.name,threshold){
  #bmp.file version
  #recommend threshold = 150~240
  #If threshold value decreases, then the number of black dots decreases.
  #file.name format = character
  #threshold format = numeric
  if(is.character(file.name)==0){
    print("Put character format in file.name")
  }
  else{
    file<-read.bmp(gsub(" ","", paste(file.name,".bmp")))
    
    arr<-array(NA,dim=c(dim(file)[1],dim(file)[2],dim(file)[3]))
    
    
    ### thresholding 
    ### black goes to more black and white goes more white.
    for(i in 1:dim(file)[1]){
      for(j in 1:dim(file)[2]){
        if(all(file[i,j,]>threshold)){arr[i,j,]=c(255,255,255)}
        else{arr[i,j,]=c(0,0,0)}
      }
    }
    
    raster.file<-as.raster(arr,max=255)
    
    x.list<-c()
    y.list<-c()
    
    
    ### #000000 is black color.
    k=1
    i=1
    e=0
    while(i <= (dim(raster.file)[1]*dim(raster.file)[2])){
      if(raster.file[i]=="#000000"){
        q=(i-1)%/%(dim(raster.file)[1])
        x=q+1
        y=(q+1)*(dim(raster.file)[1])-(i-1)
        x.list[k]=x
        y.list[k]=y
        k=k+1
      }
      else{e=e+1}
      i=i+1
    }
    
    return(list(x.list=x.list, y.list=y.list))
  }
}
setwd('/Users/cherry/Desktop')
pikachu<-dot.image("blackpika",160)
plot(pikachu$x.list,pikachu$y.list,type='p', pch=".")
length(pikachu$x.list)

