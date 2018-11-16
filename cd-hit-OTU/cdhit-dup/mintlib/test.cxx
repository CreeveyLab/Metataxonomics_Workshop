
#include<stdio.h>
#include"minString.hxx"
#include"minArray.hxx"
#include"minMap.hxx"

int main( int argc, char *argv[] )
{
	Min::String text( "The quick brown fox jumps over the lazy dog" );
	printf( "md5: %s\n", text.MD5().Data() );
	text = "";
	printf( "md5: %s\n", text.MD5().Data() );

	Min::Array<Min::String> texts;
	texts.Append( "test" );
	texts.Append( "test2" );
	printf( "%s\n", texts[1].Data() );
	return 0;
}
