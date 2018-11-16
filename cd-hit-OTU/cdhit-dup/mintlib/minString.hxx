//=================================================================
// This file is a part of the Minimum Template Library.
// By Limin Fu (phoolimin@gmail.com, lmfu@ucsd.edu)
//================================================================= 

#ifndef __MIN_STRING_HXX__
#define __MIN_STRING_HXX__

#include<string.h>
#include"minBase.hxx"

BEGIN_NS_MIN

class String
{
	char *data;
	int   size;
	int   bufsize;
	
	public:
	
	String( char ch );
	String( const char *s = NULL );
	String( const char *s, int n );
	String( const String & other );
	
	~String(){ Clear(); }
	const char* Data()const{ return data; }
	int   Size()const{ return size; }
	int   BufferSize()const{ return bufsize; }
	
	void Clear();
	void Reset();
	void Reserve( int n );
	void Resize( int n, char ch=0 );
	void SetBytes( const char *s, int n );
	void AppendBytes( const char *s, int n );
	
	void Erase( int start, int count );
	void Chop( int count = 1 );
	void ToUpper();
	void ToLower();
	
	int Find( char ch )const;
	int Find( const String & s )const;
	
	String& operator=( const char ch );
	String& operator=( const char *s );
	String& operator=( const String & s );
	
	void operator+=( const char ch );
	void operator+=( const char *s );
	void operator+=( const String &s );
	
	char& operator[]( int i ){ return data[i]; }
	char  operator[]( int i )const{ return data[i]; }
	
	String FromHex()const;
	String ToHex()const;
	String MD5()const;

	static String FromNumber( int i );
};

String operator+( const String & left, const String & right );
bool operator==( const String & left, const String & right );
bool operator!=( const String & left, const String & right );

unsigned int MakeHash( const String & s, int length=0 );
unsigned int MakeHash( const String & s, int offset, int length );

int Compare( const String & left, const String & right );
int ComparePrefix( const String & left, const String & right );
int ComparePrefix( const String & left, int lstart, const String & right, int rstart, int max=0 );

END_NS_MIN

#endif
