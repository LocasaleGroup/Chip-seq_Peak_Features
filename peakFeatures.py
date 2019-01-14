#convert from the cpp file from :https://github.com/LocasaleLab/H3K4me3_MET_Restriction/blob/master/ChIP-seq/Codes/MouseLiver/PeakFeatures.cpp
#Calculate peak height, area, width, entropy, skewness, kurtosis and 5th moment
from math import sqrt, pow, log
def calArea(values):
    return sum(values)

def calCenter(values,area):
    tmp = map(lambda x:x*values[x]/area, range(len(values)))
    return sum(list(tmp))

def calStd(values,area,center):
    tmp = map(lambda x:(values[x]/area)*(x-center)*(x-center),range(len(values)))
    return sqrt(sum(list(tmp)))

def calMoment(values,area,center,std,n):
    tmp = map(lambda x:pow((x-center)/std,n)*values[x]/area, range(len(values)))
    return sum(list(tmp))

def calEntropy(values,area,std):
    tmp = map(lambda x:std*values[x]/area, range(len(values)))
    tmp = filter(lambda x:x!=0, tmp)
    tmp = map(lambda x:x*log(x)/std,tmp)

    return -sum(list(tmp))

def calHeight(values):
    return max(values)

def peakFeaturesHeader():
    return ["height", "area", "width", "entropy", "skewness", "kurtosis", "5th_moment"]

def peakFeatures(values):

    area = calArea(values)
    if area <= 0:
        return [None]*7

    h = calHeight(values)
    center = calCenter(values,area)
    std = calStd(values,area,center)
    entropy = calEntropy(values,area,std)
    retValues = [h,area,std,entropy]
    for i in range(3,6):
        retValues.append(calMoment(values,area,center,std,i))
    return retValues
