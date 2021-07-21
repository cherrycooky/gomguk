library(EBImage)
#setwd('/Users/grrrr/Desktop/')
#charmer<-image_read("https://pngimg.com/uploads/pokemon/pokemon_PNG154.png")
#image_write(charmer,path="/Users/grrrr/Desktop/charmer.jpg",format="jpg")

image<-readImage("pika.jpg")
#image <- resize(image, w=128)
display(image)
display(image,method="raster")
str(image)
#colorMode(image) = Grayscale
#display(image)

#fHigh -> High-pass Laplacian filter 
fHigh<-matrix(1,nc=3,nr=3)

fHigh[2,2]<- -8

image.fHigh<-filter2(image,fHigh)

display(image.fHigh,method="raster")

writeImage(image.fHigh,file="dd.jpg") # jpg파일로 저장
fimage<-readJPEG("dd.jpg") # 저장한 이미지 불러오기

str(fimage) # rgb 형식으로 0~1 범위로 저장됨을 확인

dm <- dim(fimage) # 그림의 차원 할당



x=rep(1:dm[2], each=dm[1])

y=rep(dm[1]:1, dm[2])



rdata<-as.vector(fimage[,,1]*255)

gdata<-as.vector(fimage[,,2]*255)

bdata<-as.vector(fimage[,,3]*255)



bdist<-sqrt(rdata^2+gdata^2+bdata^2)

wdist<-sqrt((rdata-255)^2+(gdata-255)^2+(bdata-255)^2)

# 검정색과 하얀색까지의 거리 (유클리드)

# bdist : 검정색(0,0,0) 까지 거리

# wdist : 흰색(255,255,255)까지 거리


#thresholding 0 ~ 1
threshold=0.6
bw<-bdist<(wdist*threshold)

# bdist와 wdist 거리 비교 
#
# wdist에 원하는 상수를 곱해 준다 ex: 0.5 / 0.1 / 0.05 

result<-data.frame(cbind(x,y,rdata,gdata,bdata,bdist,wdist,bw))
head(result)

c <- result[which(result[,8]==0),]
# 테두리라고 판별된 ( 상대적으로 검은색과 거리가 가까운) 것만 저장
plot(y ~ x, data=c,pch=".") 

# extract x and y coordinates
x.list<-c()
y.list<-c()
for(i in 1:nrow(c)){
  x.list[i]=c[i,1]
}
for(i in 1:nrow(c)){
  y.list[i]=c[i,2]
}

plot(x.list,y.list,pch=20)


