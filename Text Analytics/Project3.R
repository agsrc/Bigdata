library(tm)
library(quanteda)
library(corpustools)
library(stringi)
library(collections)
library(wordnet)
library(zipfR)
library(languageR)
library(wordcloud)
library(stringr)
library(tidytext)
library(topicmodels)
library(dplyr)
library(ggplot2)

label_pos <- function(chapter){
  #Note: assumes WordNet dictionary at /usr/local/Cellar/wordnet/3.1/dict
  initDict(pathData = "/usr/local/Cellar/wordnet/3.1/dict")
  chapterClean <- gsub("[^[:alpha:][:space:]]*", "", chapter)
  chapterDTM <- DocumentTermMatrix(VCorpus(VectorSource(chapterClean)))
  chapterTerms <- chapterDTM$dimnames$Terms
  for(chapterTerm in chapterTerms){
    termLength <- nchar(chapterTerm)
    #Find all terms at least 5 characters in length
    if(termLength > 4){
      #This is annoying- synsets are grouped by POS, and in order to retrieve
      # the relevant synset from Wordnet we need to provide the POS tag.
      # So, we find the number of index terms for each possible POS type 
      # (i.e.: "ADJECTIVE", "ADVERB", "NOUN", or "VERB"), and use which ever
      # contains the most
      filter <-  wordnet::getTermFilter("ContainsFilter", chapterTerm, TRUE)
      adjectiveIndexTerms <- wordnet::getIndexTerms("ADJECTIVE", 5000, filter)
      adverbIndexTerms <- wordnet::getIndexTerms("ADVERB",5000,filter)
      nounIndexTerms <- wordnet::getIndexTerms("NOUN",5000,filter)
      verbIndexTerms <- wordnet::getIndexTerms("VERB", 5000, filter)
      
      maximalTerms <- length(adjectiveIndexTerms)
      maximalType <- "ADJECTIVE"
      indexTerms <- adjectiveIndexTerms
      if(length(adverbIndexTerms) > maximalTerms){
        maximalTerms <- length(adverbIndexTerms)
        maximalType <- "ADVERB"
        indexTerms <- adverbIndexTerms
      }
      if(length(nounIndexTerms) > maximalTerms){
        maximalTerms <- length(nounIndexTerms)
        maximalType <- "NOUN"
        indexTerms <- nounIndexTerms
      }
      if(length(verbIndexTerms) > maximalTerms){
        maximalTerms <- length(verbIndexTerms)
        maximalType <- "VERB"
        indexTerms <- verbIndexTerms
      }
      if(maximalTerms > 1 && (maximalType == "NOUN" || maximalType == "VERB")){
        synsets <- getSynsets(indexTerms[[1]])
        synset <- synsets[[1]]
        pos <- .jcall(synset, "Lcom/nexagis/jawbone/PartOfSpeech;", "getPOS")
        posString <- J("java.lang.String")$valueOf(pos)
        print(sprintf("Word: %s, POS tag: %s", chapterTerm, posString))
      }
    }
  }
}

find_n_longest_sentences <- function(jekyllVector, 
                                     n, 
                                     findLongestWord = FALSE, 
                                     findSentenceLength = FALSE,
                                     findShortestSentence = FALSE){
  #Parse sentences
  sentenceVector <- stri_split_boundaries(paste(jekyllVector, collapse=" "), type="sentence")
  #Instantiate queue for finding longest sentences in jekyllVector
  #Note that a priority queue (assuming the use of a fibonacci heap) will be much
  # faster than finding the sorting the entire list of sentences
  longestSentences <- PriorityQueue$new()
  shortestSentence <- ""
  shortestSentenceLength <- Inf
  
  #Iterate over every sentence and insert into priorty queue
  for(i in 1:length(sentenceVector[[1]])){
    sentence <- sentenceVector[[1]][i]
    #Find the number spaces in the sentence, as this corresponds to the number of words
    numSpaces <- lengths(regmatches(sentence, gregexpr(" ", sentence)))
    #Insert into the queue, where the length of the sentence is the priority
    longestSentences$push(sentence, priority = numSpaces)
    #Keep tally of the shortest sentence
    if(numSpaces < shortestSentenceLength){
      shortestSentenceLength <- numSpaces
      shortestSentence <- sentence
    }
  }
  
  #Print n longest sentences
  for(i in 1:n){
    if(findSentenceLength){
      sentence <- longestSentences$pop()
      cat(sprintf("Sentence: \"%s\", which contains %d words \n", sentence, 
            lengths(regmatches(sentence, gregexpr(" ", sentence))) + 1))
    } else{
      cat(sprintf("Sentence: \"%s\" \n",longestSentences$pop()))
    }
  }
  
  #If specified, print the shortest sentence 
  if(findShortestSentence){
    if(shortestSentenceLength > 1){
      shortestSentenceLength <- shortestSentenceLength + 1
    }
    cat(sprintf("Shortest Sentence: \"%s\", which contains %d words \n", shortestSentence, 
                shortestSentenceLength))
  }
  
  if(findLongestWord){
    jekyllDTM <-  DocumentTermMatrix(VCorpus(VectorSource(jekyllVector)))
    longestWordLength <- -99
    longestWord <- ""
    jekyllTerms <- jekyllDTM$dimnames$Terms
    for(jekyllTerm in jekyllTerms){
      termLength <- nchar(jekyllTerm)
      if(termLength > longestWordLength){
        longestWordLength <- termLength
        longestWord <- jekyllTerm
      }
    }
    print(sprintf("Longest word: %s, which is of length: %d \n", longestWord, longestWordLength))
  }
}

print_ngrams <- function(tokens, n){
  ngrams <- quanteda::tokens_ngrams(tokens, n = n, concatenator = " ")
  for(i in 1:length(ngrams)){
    for(ngram in ngrams[[i]]){
      ngramVector <- strsplit(ngram, " ")
      if(n == 2){
        if(nchar(ngramVector[[1]][1]) > 6 && nchar(ngramVector[[1]][2]) > 6){
          print(ngram)
        }
      } else{
        if(nchar(ngramVector[[1]][1]) > 6 && nchar(ngramVector[[1]][2]) > 6 && nchar(ngramVector[[1]][3]) > 6){
          print(ngram)
        }
      }
    }
  }
}

build_dendogram <- function(tokens, chapterTitle){
  #First, find the position of each word in the chapter
  #Build position matrix
  positionMatrix = matrix(nrow=length(tokens), ncol=length(tokens))
  
  for(i in 1:nrow(positionMatrix)){
    for(j in 1:ncol(positionMatrix)){
      positionMatrix[i,j] = abs(i - j)
    }
  }
  
  #Convert position matrix to object of type dist
  distMatrix <- as.dist(positionMatrix)
  
  fit <- hclust(distMatrix, method = "ward.D2")
  plot(fit, labels=tokens, main=sprintf("%s Cluster Dendrogram", chapterTitle, hang = -3))
}

build_word_cloud <- function(chapter){
  #Clean the chapter
  chapterClean <- gsub("[^[:alpha:][:space:]]*", "",chapter)
  chapterClean <- VCorpus(VectorSource(chapterClean))
  chapterClean <- tm_map(chapterClean, content_transformer(tolower))
  chapterClean <- tm_map(chapterClean, removeNumbers)
  chapterClean <- tm_map(chapterClean, removeWords, stopwords("english"))
  chapterClean <- tm_map(chapterClean, removePunctuation)
  chapterClean <- tm_map(chapterClean, stripWhitespace)
  #Build a term document matrix
  tdm <- TermDocumentMatrix(VCorpus(VectorSource(chapterClean)))
  #Determine word frequencies
  frequencies <- sort(rowSums(as.matrix(tdm)),decreasing=TRUE)
  worldcloudData <- data.frame(word = names(frequencies),freq=frequencies)
  wordcloud(words = worldcloudData$word, freq = worldcloudData$freq, min.freq = 3,
            max.words=200, random.order=FALSE, colors=brewer.pal(7, "Dark2"))
}


#Set working directory
setwd('/Users/kevintyler/Documents/Programming/csci6444')

#Store as volatile corpus
#Note: we saved, then re-opened the original file in UTF-8 encoding in a text editor
jekyllCorpus <- VCorpus(DirSource("./jekyll", mode="text",encoding ="UTF-8" ))

#Store unabridged text
fullText <- jekyllCorpus[[1]]$content

#Metadata is stored on lines 1 - 62 and 2583 - 2948, 
# so we iterate from lines 63 to 2582 to get the actual text
fullJekyllVector <- c()
chapter1JekyllVector <- c()
chapter1TokenVector <- c()
chapter2JekyllVector <- c()
chapter2TokenVector <- c()
chapter3JekyllVector <- c()
chapter3TokenVector <- c()
chapter4JekyllVector <- c()
chapter4TokenVector <- c()
chapter5JekyllVector <- c()
chapter5TokenVector <- c()
chapter6JekyllVector <- c()
chapter6TokenVector <- c()
chapter7JekyllVector <- c()
chapter7TokenVector <- c()
chapter8JekyllVector <- c()
chapter8TokenVector <- c()
chapter9JekyllVector <- c()
chapter9TokenVector <- c()
chapter10JekyllVector <- c()
chapter10TokenVector <- c()

for (i in 63:2582){
  line <- fullText[i]
  fullJekyllVector <- c(fullJekyllVector, line)
  if(i >= 62 && i <= 285){
    chapter1JekyllVector <- c(chapter1JekyllVector, line)
    splitVector <- str_split(gsub("[^[:alpha:][:space:]]*", "", line), " ")
    for(j in 1:length(splitVector[[1]])){
      chapter1TokenVector <- c(chapter1TokenVector, splitVector[[1]][j])
    }
  } else if(i >= 291 && i <= 589){
    chapter2JekyllVector <- c(chapter2JekyllVector, line)
    splitVector <- str_split(gsub("[^[:alpha:][:space:]]*", "", line), " ")
    for(j in 1:length(splitVector[[1]])){
      chapter2TokenVector <- c(chapter2TokenVector, splitVector[[1]][j])
    }
  } else if(i >= 595 && i <= 677){
    chapter3JekyllVector <- c(chapter3JekyllVector, line)
    splitVector <- str_split(gsub("[^[:alpha:][:space:]]*", "", line), " ")
    for(j in 1:length(splitVector[[1]])){
      chapter3TokenVector <- c(chapter3TokenVector, splitVector[[1]][j])
    }
  } else if(i >= 683 && i <= 838){
    chapter4JekyllVector <- c(chapter4JekyllVector, line)
    splitVector <- str_split(gsub("[^[:alpha:][:space:]]*", "", line), " ")
    for(j in 1:length(splitVector[[1]])){
      chapter4TokenVector <- c(chapter4TokenVector, splitVector[[1]][j])
    }
  } else if(i >= 844 && i <= 1026){
    chapter5JekyllVector <- c(chapter5JekyllVector, line)
    splitVector <- str_split(gsub("[^[:alpha:][:space:]]*", "", line), " ")
    for(j in 1:length(splitVector[[1]])){
      chapter5TokenVector <- c(chapter5TokenVector, splitVector[[1]][j])
    }
  } else if(i >= 1032 && i <= 1165){
    chapter6JekyllVector <- c(chapter6JekyllVector, line)
    splitVector <- str_split(gsub("[^[:alpha:][:space:]]*", "", line), " ")
    for(j in 1:length(splitVector[[1]])){
      chapter6TokenVector <- c(chapter6TokenVector, splitVector[[1]][j])
    }
  } else if(i >= 1171 && i <= 1230){
    chapter7JekyllVector <- c(chapter7JekyllVector, line)
    splitVector <- str_split(gsub("[^[:alpha:][:space:]]*", "", line), " ")
    for(j in 1:length(splitVector[[1]])){
      chapter7TokenVector <- c(chapter7TokenVector, splitVector[[1]][j])
    }
  } else if(i >= 1236 && i <= 1716){
    chapter8JekyllVector <- c(chapter8JekyllVector, line)
    splitVector <- str_split(gsub("[^[:alpha:][:space:]]*", "", line), " ")
    for(j in 1:length(splitVector[[1]])){
      chapter8TokenVector <- c(chapter8TokenVector, splitVector[[1]][j])
    }
  } else if(i >= 1722 && i <= 1989){
    chapter9JekyllVector <- c(chapter9JekyllVector, line)
    splitVector <- str_split(gsub("[^[:alpha:][:space:]]*", "", line), " ")
    for(j in 1:length(splitVector[[1]])){
      chapter9TokenVector <- c(chapter9TokenVector, splitVector[[1]][j])
    }
  } else {
    chapter10JekyllVector <- c(chapter10JekyllVector, line)
    splitVector <- str_split(gsub("[^[:alpha:][:space:]]*", "", line), " ")
    for(j in 1:length(splitVector[[1]])){
      chapter10TokenVector <- c(chapter10TokenVector, splitVector[[1]][j])
    }
  }
}

#Print the the longest sentences in the entire corpus
print("Printing the top ten longest sentences in the entire document...")
find_n_longest_sentences(fullJekyllVector,10)

print("Printing the longest word, longest sentence, and shortest sentence in each chapter...")
#Find the longest sentence and longest word for each sentence
find_n_longest_sentences(chapter1JekyllVector,1, TRUE, TRUE, TRUE)
find_n_longest_sentences(chapter2JekyllVector,1, TRUE, TRUE, TRUE)
find_n_longest_sentences(chapter3JekyllVector,1, TRUE, TRUE, TRUE)
find_n_longest_sentences(chapter4JekyllVector,1, TRUE, TRUE, TRUE)
find_n_longest_sentences(chapter5JekyllVector,1, TRUE, TRUE, TRUE)
find_n_longest_sentences(chapter6JekyllVector,1, TRUE, TRUE, TRUE)
find_n_longest_sentences(chapter7JekyllVector,1, TRUE, TRUE, TRUE)
find_n_longest_sentences(chapter8JekyllVector,1, TRUE, TRUE, TRUE)
find_n_longest_sentences(chapter9JekyllVector,1, TRUE, TRUE, TRUE)
find_n_longest_sentences(chapter10JekyllVector,1,TRUE, TRUE, TRUE)

#Use WordNet to mark the parts of speech for the first chapter for nouns and verbs having a length of 5 or greater.
label_pos(chapter1JekyllVector)

#Analyze word frequency using functions from package zipfR.
#Note- this follows the example provided by Stefan Ever in
# documentation of zipfR
jekyllSpc <- text2spc.fnc(fullJekyllVector)
zm <- lnre("zm",jekyllSpc)
zm.spc <- lnre.spc(zm,N(jekyllSpc))
plot(jekyllSpc,zm.spc, title="Word Frequency: Dr. Jekyll and Mr. Hyde")
EV(zm,1e+8)
EVm(zm,1,1e+8)
ext.vgc <- lnre.vgc(zm,(1:100)*1e+5,m.max=1)
plot(ext.vgc,N0=N(zm),add.m=1, main = "Dr. Jekyll and Mr. Hyde Expected V and V1 growth curves")

#Generate bigrams and trigrams for all words whose length is greater than 6 characters in Chapter 1.
chapter1Clean <- gsub("[^[:alpha:][:space:]]*", "", chapter1JekyllVector)
chapter1Tokens <- quanteda::tokens(corpus(chapter1Clean), removePunct = TRUE, removeNumbers = TRUE)
#Print bigrams
print_ngrams(chapter1Tokens,2)
#Print trigrams
print_ngrams(chapter1Tokens,3)

#Build dendrograms
build_dendogram(chapter1TokenVector, "Chapter 1")
build_dendogram(chapter2TokenVector, "Chapter 2")
build_dendogram(chapter3TokenVector, "Chapter 3")
build_dendogram(chapter4TokenVector, "Chapter 4")
build_dendogram(chapter5TokenVector, "Chapter 5")
build_dendogram(chapter6TokenVector, "Chapter 6")
build_dendogram(chapter7TokenVector, "Chapter 7")
build_dendogram(chapter8TokenVector, "Chapter 8")
build_dendogram(chapter9TokenVector, "Chapter 9")
build_dendogram(chapter10TokenVector, "Chapter 10")

#Build wordclouds
build_word_cloud(chapter1JekyllVector)
build_word_cloud(chapter2JekyllVector)
build_word_cloud(chapter3JekyllVector)
build_word_cloud(chapter4JekyllVector)
build_word_cloud(chapter5JekyllVector)
build_word_cloud(chapter6JekyllVector)
build_word_cloud(chapter7JekyllVector)
build_word_cloud(chapter8JekyllVector)
build_word_cloud(chapter9JekyllVector)
build_word_cloud(chapter10JekyllVector)

#Miscellaneous corpustools operations
chapter1TC <- corpustools::create_tcorpus(chapter1JekyllVector)
chapter1TC$get()
chapter1TC$subset(subset = docfreq_filter(token, min=2))

#Miscellaneous stringi operations
chapter1UniqueStrings <- stri_unique(chapter1JekyllVector)
chapter1EmptyStringRemoved <- stri_remove_empty(chapter1JekyllVector)

#Miscellaneous tidytext operations
chapter1Doc <- data_frame(document = chapter1JekyllVector)
chapter1Doc %>% unnest_tokens(word, document)
#Instantiate the tidy class
chapter1Tidy <- tidy(chapter1JekyllVector)
#Clean the data
chapterClean <- gsub("[^[:alpha:][:space:]]*", "", chapter1JekyllVector)
#Remove all zero topics
chapterDTM <- DocumentTermMatrix(VCorpus(VectorSource(chapterClean)))
chapterLDA <- LDA(chapterDTM, 5)
td_lda <- tidy(chapterLDA)
td_lda_filtered <- td_lda %>%
  filter(beta > .004) %>%
  mutate(term = reorder(term, beta))
#Plot the top terms in each topic of chapter 1
ldaPlot <- ggplot(td_lda_filtered, aes(term, beta)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ topic, scales = "free") +
  theme(axis.text.x = element_text(angle = 90, size = 15))
#Add a title
ldaPlot + ggtitle("Dr. Jekyll Chapter 1 Topic Analysis", 
                  theme(plot.title = element_text(lineheight=.8, face="bold")))