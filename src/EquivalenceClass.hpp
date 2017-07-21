/*

Copyright (c) 2013  Ghassen Hamrouni

Abstract:

This module provides models an equivalence class using
disjoint set forest.

Author:

Ghassen Hamrouni <ghamrouni.iptech@gmail.com> 24-07-2013

Revision History:

*/

#pragma once

#ifndef ___CORE_DISJOINT_SET_FORESTS_H__
#define ___CORE_DISJOINT_SET_FORESTS_H__

#include <memory>

template<class T>
class EquivalenceClass
{
public:

	T data;
	EquivalenceClass<T> *parent;
	int rank;

	explicit EquivalenceClass(T data) {
		this->data = data;
		rank = 0;
		parent = this;
	}

	virtual ~EquivalenceClass() { }

	EquivalenceClass(const EquivalenceClass<T>& p)
	{
		this->data = p.data;
		rank = p.rank;
		parent = p.parent;
	}

	EquivalenceClass<T>& operator=(const EquivalenceClass<T>& p)
	{
		this->data = p.data;
		rank = p.rank;
		parent = p.parent;
	}

	EquivalenceClass<T>* findRepresentative() {
		if (parent == this) 
			return this;

		EquivalenceClass<T>* rep = parent->findRepresentative();

		// Path compression
		parent = rep;

		return rep;
	}

	// union by rank
	void merge(EquivalenceClass<T>* e1) {
		if (e1->rank > rank)
		{
			parent = e1;
		} else { 

			e1->parent = this;

			if (e1->rank == rank)
			{
				rank += 1;
			}
		}
	}

	T get() const {
		return data;
	}

	void set(const T& data){
		this.data = data;
	}

	bool operator==(EquivalenceClass<T> e1)
	{
		EquivalenceClass<T>* rep1 = e1.findRepresentative();
		EquivalenceClass<T>* rep2 = findRepresentative();

		return (rep1->parent == rep2->parent);
	}
};

#endif // ___CORE_DISJOINT_SET_FORESTS_H__