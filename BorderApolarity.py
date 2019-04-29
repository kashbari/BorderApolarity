# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 20:16:51 2019

@author: kashb
"""


def part(n, k):
	def _part(n, k, pre):
		if n <= 0:
			return []
		if k == 1:
			if n <= pre:
				return [[n]]
			return []
		ret = []
		for i in range(min(pre, n), 0, -1):
			ret += [[i] + sub for sub in _part(n-i, k-1, i)]
		return ret
	return _part(n, k, n)

def part1(n,k):
	p = []
	for i in range(1,k+1):
		p += part(n,i)
	return p


