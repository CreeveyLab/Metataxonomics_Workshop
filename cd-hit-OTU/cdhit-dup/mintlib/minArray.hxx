//=================================================================
// This file is a part of the Minimum Template Library.
// By Limin Fu (phoolimin@gmail.com, lmfu@ucsd.edu)
//================================================================= 

#ifndef __MIN_ARRAY_HXX__
#define __MIN_ARRAY_HXX__

#include <string.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>
#include <new>
#include "minBase.hxx"

BEGIN_NS_MIN

template<typename item_t>
class Array
{
	item_t *items;
	item_t *buf;
	int     size;
	int     bufsize;
	
public:
	Array( int n=0, item_t it=item_t() ){
		items = buf = NULL;
		size = bufsize = 0;
		Resize( n, it );
	}
	Array( const Array & other ){
		items = buf = NULL;
		size = bufsize = 0;
		this->operator=( other );
	}
	~Array(){
		Clear();
	}
	
	int Size()const{ return size; }
	item_t& operator[]( int i ){ return items[i]; }
	item_t  operator[]( int i )const{ return items[i]; }
	
	Array& operator=( const Array & other ){
		int i;
		Clear();
		size = other.size;
		bufsize = other.bufsize;
		items = buf = (item_t*) calloc( bufsize, sizeof(item_t) );
		for(i=0; i<size; i++) new (items+i) item_t( other.items[i] );
		return *this;
	}
	
	void Swap( Array & other ){
		item_t *it = items, *bu = buf;
		int s = size, bs = bufsize;
		items = other.items;
		buf = other.buf;
		size = other.size;
		bufsize = other.bufsize;
		other.items = it;
		other.buf = bu;
		other.size = s;
		other.bufsize = bs;
	}
	void Clear(){
		int i;
		for(i=0; i<size; i++) items[i].~item_t();
		if( buf ) free( buf );
		items = buf = NULL;
		size = bufsize = 0;
	}
	void ResetBuffer(){ // remove buffer from front:
		if( items == buf ) return;
		memmove( buf, items, size*sizeof(item_t) );
		items = buf;
	}
	void Reserve( int n, bool extra=false ){
		int front = (intptr_t) (items - buf);
		if( bufsize >= (n + front) ) return; // enough space at back
		if( bufsize >= n ){ // front > 0
			ResetBuffer();
			return;
		}
		void *old = buf;
		if( extra ) n += n/5;
		buf = (item_t*) malloc( n*sizeof(item_t) );
		memmove( buf, items, size*sizeof(item_t) );
		items = buf;
		bufsize = n;
		free( old );
	}
	void Resize( int n, item_t it=item_t() ){
		int i;
		ResetBuffer();
		if( bufsize != n ){
			items = buf = (item_t*) realloc( buf, n*sizeof(item_t) );
			bufsize = n;
		}
		for(i=size; i<n; i++) new (items+i) item_t( it );
		size = bufsize = n;
	}
	item_t& Front(){
		assert( size );
		return items[0];
	}
	item_t& Back(){
		assert( size );
		return items[size-1];
	}
	void PushFront( const item_t & it ){
		int front = (intptr_t) (items - buf);
		if( front ){
			items -= 1;
		}else{
			front = bufsize/5 + 5;
			bufsize += front;
			buf = (item_t*) malloc( bufsize*sizeof(item_t) );
			memmove( buf + front, items, size*sizeof(item_t) );
			free( items );
			items = buf + front - 1;
		}
		new (items) item_t( it );
		size += 1;
	}
	void PushBack( const item_t & it ){
		Reserve( size + 1, true );
		new (items+size) item_t( it );
		size += 1;
	}
	void Append( const item_t & it ){ PushBack( it ); }
	void Erase( int start, int n=-1 ){
		if( size == 0 ) return;
		if( start < 0 or start >= size ) return;
		if( n < 0 ) n = size - start;
		n += start;
		if( n > size ) n = size;
		
		int i;
		for(i=start; i<n; i++){
			items[i].~item_t();
			memset( items+i, 0, sizeof(item_t) );
		}
		for(i=0; i<size-n; i++){
			memcpy( items+start+i, items+n+i, sizeof(item_t) );
		}
		size -= n - start;
	}
	void PopBack(){ if( size ) Erase( size-1, 1 ); }
};

END_NS_MIN

#endif
