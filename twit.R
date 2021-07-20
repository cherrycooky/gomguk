library(base64enc)
library(rtweet)
library(twitteR)

#tweeteR
consumer_key   <- "7JuRFQqpzCN0BBBkOur7DukW5"
consumer_secret<- "OqiuAEimHUU2DaVgw2RJEiP1Tt78Y6RHiZqxxToMroM7980X8r"
access_token   <- "1417360319932862469-ngL3GV3iFRv4pfmhAodmTRCn5gWmQP"
access_secret  <- "OTBGQttEy96QdKW5QIz4fAQKgD2Fp09eJC7714Wv3SKkf"

twitter_token <- create_token(app = "gomgukwow", consumer_key = consumer_key, consumer_secret = consumer_secret, access_token = access_token, access_secret = access_secret)


post_tweet(status="fourth in a thread", token = twitter_token)

