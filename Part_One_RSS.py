##Python Assignment 3 Part 1
##Sarah McCartan
##15206644
##COMP41680

#The goal is to construct a corpus of documents with 2
#classes, by collecting data from 2 different RSS news feeds. See the Appendix of this
#assignment for more details on RSS feeds.
#Tasks to be completed:
#1. Identify 2 different online news RSS feeds. It is recommended that these RSS
#feeds are on distinct topics (e.g. business and sport; movies and music etc), and
#that the items in the feed include descriptions with snippets of article text.
#2. For each feed (i.e. class):
#(i) Collect a reasonable number of items (i.e. articles). This may require rerunning
#the collection process over a period of time.
#(ii) Parse the news feeds, and extract the text of all item description. These
#texts will represent the documents for this class.
#(iii) Write the documents to a file format of your choice.

#import all the libraries necessry to run this script
import urllib.request
import csv
import pandas as pd
import xml.etree.ElementTree as tree

#First RSS feed identified is theguardian.com/fashion/rss

rss_fashion=[] #creating empty list for the Fashion RSS Feed

url = "http://www.theguardian.com/fashion/rss"
response = urllib.request.urlopen(url) #action open url 
data=tree.parse(response) #parse xml data
    
for entry in data.iterfind("channel/item"): # iterate through xml tag channel sub-tag item to find relevant data 
    description = entry.findtext("description")#pull out the description text
    article_id = entry.findtext("guid") #pull out unique id for each article
    rss_fashion.append(description) # append to list
    found = False #creating boolean operation for ensuring addition of only unique articles
    for article_id in rss_fashion: 
        if article_id == article_id:
            found = True
            break
        if not found: #if article id is not found then it is unique and needs to be added to the list
            rss_fashion.append(description)

fashion_label=["fashion"]*len(rss_fashion) #create list of fashion class labels
fashion=pd.DataFrame(rss_fashion) #convert list to dataframe
fashiondf=fashion.assign(class_label = fashion_label) #add to dataframe
#fashion.columns["desc", "class_label"] #assign column names
fashiondf.to_csv("fashion-raw-rss.csv") #convert df to csv file

#Second RSS feed identified is guardian.com/sports/rss

rss_sport=[] #creating empty list for the Sport RSS Feed

url = "http://www.theguardian.com/sport/rss"
response = urllib.request.urlopen(url) #action open url 
data=tree.parse(response) #parse xml data
    
for entry in data.iterfind("channel/item"): # iterate through xml tags channel sub-tag item find data 
    description = entry.findtext("description")#description text
    article_id = entry.findtext("guid") #unique id for each article
    rss_sport.append(description) # append to list
    found = False #creating boolean operation for ensuring addition of only new articles
    for article_id in rss_sport: 
        if article_id == article_id:
            found = True
            break
        if not found: #if article id is not found then it is unique and needs to be added to the list
            rss_sport.append(description)
			
sport_label=["sport"]*len(rss_sport)
sport=pd.DataFrame(rss_sport) #convert list to Pandas DataFrame
sportdf=sport.assign(class_label=sport_label)
#sport.columns["desc", "class_label"]
sportdf.to_csv("sports-raw-rss.csv") #convert series to csv file