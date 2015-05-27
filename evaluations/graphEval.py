#!/usr/bin/env python
from __future__ import print_function,division
import urllib2, json, argparse, sys
from collections import OrderedDict

def getAlleles(url):
	alleleDict=OrderedDict()
	success=True
	alleleNumber=-1
	while success==True:
		alleleNumber+=1
		try:
			text=urllib2.urlopen(url+'/v0.6.a/alleles/'+str(alleleNumber)).read()
			json_=json.loads(text)
			id_=json_['id']
			variantSetID=json_['variantSetId']
			path=json_['path']
			segments=path['segments']
			allelePathItemList=[]
			for segment in segments:
				allelePathItem={}
				length=int(segment['length'])
				start=segment['start']
				strand=start['strand']
				base=start['base']
				pos=int(base['position'])
				seqId=int(base['sequenceId'])
				allelePathItem['seq']=seqId
				allelePathItem['strand']=strand
				allelePathItem['pos']=pos
				allelePathItem['length']=length
				allelePathItemList.append(allelePathItem)
			alleleDict[id_]=allelePathItemList
		except urllib2.HTTPError:
			success=False
	return alleleDict

def getRefOverlap(allelePathItemList,refPathList):
	overlapLength=0
	totalLength=0
	for pathItem in allelePathItemList:
		totalLength+=pathItem['length']
	for pathItem in allelePathItemList:
		segSeq=pathItem['seq']
		segStrand=pathItem['strand']
		segStart=pathItem['pos']
		segLength=pathItem['length']
		segEnd=segStart+segLength-1
		if segStrand=='NEG_STRAND':
			segEnd=segStart-segLength+1
		segStart,segEnd=sorted([segStart,segEnd])
		for refPathItem in refPathList:
			refSeq=refPathItem['seq']
			refStrand=refPathItem['strand']
			refStart=refPathItem['pos']
			refLength=refPathItem['length']
			refEnd=refStart+refLength-1
			if refStrand=='NEG_STRAND':
				refEnd=refStart-refLength+1
			refStart,refEnd=sorted([refStart,refEnd])
			if segSeq==refSeq:
				if not ((segStart<refStart and segEnd<refEnd) or (segStart>refStart and segEnd>refEnd)):
					if segStart>=refStart and segEnd<=refEnd:
						overlapLength+=segLength
					elif segStart<=refStart and segEnd>=refEnd:
						overlapLength+=refLength
					elif segStart<=refStart and segEnd>=refStart:
						overlapLength+=segEnd-refStart+1
					elif segStart<=refEnd and segEnd>=refEnd:
						overlapLength+=segEnd-refEnd+1
	refOverlapFraction=overlapLength/totalLength
	return refOverlapFraction


def parseArgs():
	parser = argparse.ArgumentParser(description="""Performs the specified evaluation on a specified graph server.
		Requires a url to be supplied from the user.""")
	parser.add_argument('url',type=str,help="""A string containing the url of the graph server to be evaluated.""")
	parser.add_argument('--align2ref', dest='eval', action='store_const', const=getRefOverlap, default=getRefOverlap,
		help="""Specifies the evaluation to be performed.  Currently only one option, so this argument is unnecessary.""")
	parser.add_argument('--ref',default='0',help="""Specifies the id of the path (a nonnegative integer) that 
		corresponds to the reference path.""")
	args = parser.parse_args()
	return args
	


if __name__ == '__main__':
	args=parseArgs()
	alleleDict=getAlleles(args.url)
	refAllele=alleleDict[args.ref]
	for id_,allele in alleleDict.items():
		refOverlap=getRefOverlap(allele,refAllele)
		print("ID: ",id_)
		print("% overlap: ",refOverlap,"\n")