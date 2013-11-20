
func E2880(i)
/* DOCUMENT j = E2880(i)
     returns the first multiple of 2880 greater or equal than i
   SEE ALSO:
 */
{
  return (1l+(long(i)-1l)/2880l)*2880l;
}


func anything2str(xvalue)
/* DOCUMENT str = anything2str(x)
     Converts an array of numbers (floats, double or long) into string array.
     
   SEE ALSO:
 */
{
  tx=typeof(xvalue);
  if(tx=="string") return xvalue;
  value = [];
  if(tx=="char"|tx=="int"|tx=="long"|tx=="short") {
    value = swrite(format="%d",long(xvalue));
  }
  if(tx=="float"|tx=="double") {
    value = swrite(format="%g",double(xvalue));
  }
  return value;
}



func readfits( fichier, err= )
/* DOCUMENT  a = readfits( fichier )
         or  a = readfits( fichier, err=1 )
The function reads a FITS file.

When using "a = readfits( fichier )", a problem with the file 'fichier' causes
a runtime error.

When using "a = readfits( fichier , err=1 )", a problem with the file
'fichier' returns a nil ([]) parameter.


*/
{
  if( err==1 ) {
    fic = open( fichier, "rb", 1);
    if( is_void(fic) ) {
      write,format="File %s was not found by function readfits.\n",fichier;
      return [] ;
    }
  }
  else {
    fic = open( fichier, "rb");
  }

  sun_primitives,fic;

  // Reading data at the begining of file (i.e. at ptr=0)
  ptr = 0;
  dataTotalByteSize = extend = 0;
  data = core_readfits( fic, ptr, dataTotalByteSize, extend );   //  dataTotalByteSize and extend are outputs

  // if extended FITS detected, then reading extensions one by one
  if( extend ) {
    data = array(pointer, extend);
    for(i=1; i<=extend; i++) {
      // setting pointer to begining of next extension
      ptr = E2880(dataTotalByteSize);  // first multiple of 2880 greater (or equal) than i
      // then reading data
      data(i) = &core_readfits( fic, ptr, dataTotalByteSize, extend );   //  dataTotalByteSize is output, extend is input
    }
  }

  return data;
}





func readFitsField(fic, ptr, &i, keyword, format)
{
  n=0;
  if(format=="%d") val = 0;
  else if(format=="%f") val=0.0;
  else if(format=="%s") val=" ";
  else val=0;
  tmp = array(char,80);
  /*
    La c est une astuce pour lire les fichiers ali .. il faut mettre
    35e4 au lieu de 35 pour que la boucle aille chercher le keyword
    super loin hors du header et des data ... je sais plus pourquoi
    .. ah si c'est parce que ali ecrit des data a la suite de son 1er
    header, et pas moi
    Mais du coup si on met le 35e4 il faut changer un truc ailleurs ?
    keske j'ai change ailleurs ? je sais plus .. merde.
    Bon tant pis je remets le 35 en place.
  */
  // while ( i<35e4 && n!=1 ) { 
  while ( i<35 && n!=1 ) { 
    _read,fic,ptr + i*80,tmp;
    i++;
    str = string(&tmp);
    n = sread(str,format=keyword+" = "+format+" ",val);
  }
  return val;
}








func core_readfits( fic, ptr, &dataTotalByteSize, &extend )
/* DOCUMENT
   <fic> file descriptor
   <ptr> number of bytes, indicating where the reading is to begin in the file
   <dataTotalByteSize> output variable telling where in the file the current data reading
                       has finished (in bytes)
   <extend> intput/output variable. When set to 0, the beginning of a stadard FITS file is
            expected to be read, and the EXTEND and NEXTEND keywords will be seached for.
            If found, then <extend> contains the number of extensions on output.
            In input, <extend> set to a non-zero value means that a FITS extension is
            expected to be read

*/
{
  // required so that values are initialized when procedure exits on error
  dataTotalByteSize = 0;

  // variable telling whether data should be read or not (they shouldnt be
  // when the parsed header describes file extensions, instead of data)
  noread = 0;
  
  // Searching for 'SIMPLE = T' keyword ...
  if( extend==0 ) {
    i = 0;
    T = readFitsField(fic, ptr, i, "SIMPLE", "%s");
    if(i!=1) {
      write,"ERROR : Can't read SIMPLE=T keyword. File is not written in FITS ";
      return [];
    }
    if( T!="T" ) {
      if(T==[]) T="[]";
      write,format="WARNING : File is not written in true FITS.\nSIMPLE = %s\nTrying to continue anyway.\n",T;
    }
  

    // Searching a possible 'EXTEND = T' keyword ... and get the number of extensions when found.
    i=1;
    if( readFitsField(fic, ptr, i, "EXTEND", "%s")=="T" ) {
      extend = readFitsField(fic, ptr, i, "NEXTEND", "%d");
      // if extensions are *really* present ...
      if( extend!=0 )
        noread = 1; // this header is a global header defining extensions: do not attempt to read any data, and skip anything related to data reading
    }
  } else {
    // Reading the extension of a FITS file
    i = 0;
    T = readFitsField(fic, ptr, i, "XTENSION", "%s");
  }
  

  // Going back to start of file, searching for BITPIX, NAXIS, etc ....
  if( noread==0 ) {
    i=1;
    bitpix = readFitsField(fic, ptr, i, "BITPIX", "%d");

    naxis = readFitsField(fic, ptr, i, "NAXIS", "%d");
    // creation of dimension_list of data array
    dimsofNaxis = array(naxis,naxis+1);
    // reading all dimensions
    dataTotalByteSize = 1;
    for(k=1; k<=naxis; k++) {
      naxKeyw = swrite(format="NAXIS%d",k);
      dataTotalByteSize *= (dimsofNaxis(k+1) = readFitsField(fic, ptr, i, naxKeyw, "%d"));
    }
    dataTotalByteSize *= abs(long(bitpix)/8l);
  
    // data type
    intype = 0;
    if( abs(bitpix)==16 ) {
      datatype = short(0);
      if( bitpix==-16 ) intype=1;   // special case for Ali ....
    }
    else if (bitpix==-32 )
      datatype = float(0);
    else if (bitpix==64 )
      datatype = long(0);
    else if (bitpix==32 )
      datatype = long(0);
    else if (bitpix==-64 )
      datatype = double(0);
    else if (bitpix==8 )
      datatype = char(0);
    else {
      write,format="ERROR in formatread.i : unknown data type BITPIX=%d \n",bitpix;
      return [];
    }

    // allocating memory for data array, with proper dimensions
    data = array(datatype, dimsofNaxis);
  }
  
  // reperage de la position du END a la fin du header
  keywd=" ";
  tmp=array(char,80);
  while ( keywd!="END" && i<1000 ) {
    _read,fic,ptr+i*80,tmp; i++; str = string(&tmp); 
    n = sread(str,format="%s ",keywd);
  }
  
  /* comptage du nombre de blocks ds le header */  
  n = 0;
  while( i>0 ) {
    i = i-36;
    n++;
  }


  if( noread==0 ) 
    _read,fic,ptr+2880*n,data;
  

  dataTotalByteSize += ptr + 2880*n;

  // special conversion for format BITPIX=-16
  if( intype==1 ) {
    data = long(data);
    nn = where(data<0);
    if( is_array(nn) )
      data(nn) += 65536;
  }
  
  return( data );
}






/*

  
__        ______  ___ _____ _____ _____ ___ _____ ____  
\ \      / /  _ \|_ _|_   _| ____|  ___|_ _|_   _/ ___| 
 \ \ /\ / /| |_) || |  | | |  _| | |_   | |  | | \___ \ 
  \ V  V / |  _ < | |  | | | |___|  _|  | |  | |  ___) |
   \_/\_/  |_| \_\___| |_| |_____|_|   |___| |_| |____/ 
                                                        
  
 */








func writefits(fichier, data, keyword, xvalue, comment=)
/* DOCUMENT writefits, fichier, data, keyword, value, comment=
   Writes <data> in the the file called <fichier>.
   Optionally, one can provide it with some keywords with the associated
   values, and comments.

   Examples:
   writefits, "myfile.fits", data;
   writefits, "myfile.fits", data, ["DATE","OBSERVER"],[fitsDate(),"Bob"], comment="just a joke"
   writefits, "myfile.fits", data, "DATE",fitsDate(), comment=["just a joke","really !"]
   
   SEE ALSO:
 */
{
  fic = open( fichier, "wb+");
  if( !fic ) {
    write,format="Can't open file %s. Exiting now. Nothing written.\n",fichier;
    error;
  }
  sun_primitives,fic;

  ptr = 0;
  if( typeof(data)=="pointer" ) {
    nextend = numberof(data);
    core_writefits, fic, ptr, nextend, [], keyword, xvalue, comment=comment;
    for(i=1;i<=nextend;i++) {
      //      write,format="extension %d at ptr=%d \n",i,ptr;
      core_writefits, fic, ptr, nextend, *data(i);
    }
  } else {
    core_writefits, fic, ptr, 0, data, keyword, xvalue, comment=comment;
  }
  close, fic;
  return 0;
}



func writeFitsField(fic, ptr, keyword, xvalue, comment)
{
  if( keyword=="COMMENT" ) {
    tmp = fill80( "COMMENT " + comment );
  } else if( keyword=="END" ) {
    tmp = fill80("END");
  } else if( keyword=="SIMPLE" || keyword=="EXTEND" ) {
    tmp = fill80( fill80(keyword,8) + "=" + pref80("T",21) + " / "+comment );    
  } else if( keyword==" " ) {
    tmp = fill80(" ");
  } else {
    value = anything2str(xvalue);
    if( is_a_number(xvalue) ) {
      // we need to write a number
      tmp = fill80( fill80(keyword,8) + "=" + pref80( value,21) + " / "+comment );
    } else {
      // the string is really a string
      tmp = fill80( fill80(keyword,8) + "= '" + value + "'" );
    }
  }
  _write, fic, ptr, strchar(tmp);
  ptr += 80;
  return ptr;
}




func core_writefits(fichier, &ptr, nextend, data, keyword, xvalue, comment=)
/* DOCUMENT writefits, fichier, data, keyword, value, comment=
   Writes <data> in the the file called <fichier>.
   Optionally, one can provide it with some keywords with the associated
   values, and comments.

   Examples:
   writefits, "myfile.fits", data;
   writefits, "myfile.fits", data, ["DATE","OBSERVER"],[fitsDate(),"Bob"], comment="just a joke"
   writefits, "myfile.fits", data, "DATE",fitsDate(), comment=["just a joke","really !"]
   
   SEE ALSO:
 */
{
  bitpix=0;
  value = anything2str( xvalue );  // converts xvalue into a string

  // if data==[], no data should be written: only the header of extensions should be written
  if( nextend==0 ) {
    status="simple+data";
    if( data==[] )
      return -1;
  }
  else {
    if( data==[] )
      status="global-extension";
    else
    status="extension+data";
  }

  if( status!="global-extension" ) {
    // we inted to write the data ...
    datatype = typeof(data);
    if( datatype=="double" ) bitpix=-64;
    if( datatype=="float" ) bitpix=-32;
    if( datatype=="long" ) bitpix=32;
    if( datatype=="int" ) bitpix=32;
    if( datatype=="short" ) bitpix=16;
    if( datatype=="char" ) bitpix=8;
    if( !bitpix ) {
      write,format="Error : unknown data type %s \n",datatype;
      error;
    }
    dim = dimsof(data);
  }
  
  if( status=="simple+data" ) {
    ptr = writeFitsField(fic, ptr, "SIMPLE", "T", "Std FITS");
  }

  if( status=="extension+data" ) {
    ptr = writeFitsField(fic, ptr, "XTENSION", "IMAGE", "extension ...");
  }
  
  if( status!="global-extension" ) {
    ptr = writeFitsField(fic, ptr, "BITPIX", bitpix, "Yorick "+typeof(data));
    naxis = dim(1);
    ptr = writeFitsField(fic, ptr, "NAXIS", naxis, "Number of axes");
    for(i=1; i<=naxis;i++) {
      naxKeyw = swrite(format="NAXIS%d",i);
      ptr = writeFitsField(fic, ptr, naxKeyw, dim(i+1), swrite(format="Axis %d",i) );
    }
  }
  
  if( status=="global-extension" ) {
    ptr = writeFitsField(fic, ptr, "SIMPLE", "T", "Std FITS");
    ptr = writeFitsField(fic, ptr, "BITPIX", 8, "Useless");
    ptr = writeFitsField(fic, ptr, "NAXIS", 0, "No data in primary header");
    ptr = writeFitsField(fic, ptr, "EXTEND", "T", "Extended FITS");
    ptr = writeFitsField(fic, ptr, "NEXTEND", nextend, "Nber of extensions");
  }
  
  // WRITING KEYWORDS
  nk = numberof(keyword);
  nval = numberof(value);
  if( nval!=nk ) {
    write,format="WARNING : numberof keyword (%d) different from number of values (%d)\n", nk, nval;
    nk = min(nk,nval);
    write,format="WARNING : only (%d) values will be written\n", nk;
  }
  for(i=1; i<=nk; i++) {
    ptr =  writeFitsField(fic, ptr, keyword(i), value(i), " ");
  }

  // WRITING COMMENTS
  if( !is_void(comment) ) {
    nc = numberof(comment);
    for(i=1;i<=nc;i++) {
      ptr =  writeFitsField(fic, ptr, "COMMENT", 0 , comment(i));
    }
  }
  ptr =  writeFitsField(fic, ptr, "END", 0,0);

  // FILLING THE HEADER WITH ' '  ..........    // 'LE' bug du mercredi 24 Novembre 2010
  npad = E2880(ptr) - ptr;
  if( npad ) {
    tmp = array(' ', npad);
    _write, fic, ptr, tmp;
    ptr += npad;
  }

  // WRITING DATA
  if( status!="global-extension" ) {
    _write,fic,ptr,data;
    // ptr += sizeof(data);   // cree un bug sur une architecture 64 bits car sizeof(long(0))=8
    ptr += numberof(data) * abs(bitpix)/8;
  }

  
  // zero-padding ....
  npad = E2880(ptr) - ptr;
  if( npad ) {
    tmp = array(char(0), npad);
    _write, fic, ptr, tmp;
    ptr += npad;
  }
  // write,format="ptr in core_writefits = %d \n",ptr;
  return 0;
}










func retspace(i)
/* DOCUMENT str = retspace(i)
     Returns a string of <i> whitespaces.
   SEE ALSO:
 */
{
  if(i<1) return "";
  return string(&array(' ',i));
}


func fill80(str, k)
/* DOCUMENT fill80(str, k)
     Complement right part of string <str> with
     whitespaces up to a total of <k> characters.
     <k> can be omitted, default value is 80.
   SEE ALSO:
 */
{
  if( is_void(k) ) k=80;
  l = strlen(str);
  if( l>k ) {
    return strpart(str,1:k);
  }
  return str+retspace(k-l);
}


func pref80(str, k)
/* DOCUMENT fill80(str, k)
     Complement left part of string <str> with
     whitespaces up to a total of <k> characters.
     <k> can be omitted, default value is 80.
   SEE ALSO:
 */
{
  if( is_void(k) ) k=80;
  l = strlen(str);
  if( l>k ) {
    return strpart(str,l-k+1:l);
  }
  return retspace(k-l) + str;
}




func putBinaryRight(str, tmp, k0, k1)
{
  tmp( k0:k1 ) = ' ';
  l = strlen(str);
  if(l==0) return tmp;
  ka = k1-l+1;
  strc = strchar(str)(:-1);
  k0 = max(ka,k0);
  tmp( k0:k1 ) = strc(1:(k1-k0+1));
  return tmp;
}


func putBinaryLeft(str, tmp, k0, k1)
{
  tmp( k0:k1 ) = ' ';
  l = strlen(str);
  if(l==0) return tmp;
  kb = k0+l-1;
  strc = strchar(str)(:-1);
  k1 = min(kb,k1);
  tmp( k0:k1 ) = strc(1:(k1-k0+1));
  return tmp;
}





func is_a_number( xvalue )
/* DOCUMENT  is_a_number(xvalue)

     Argument <xvalue> may either be a number or a string.

     The function returns 1 when xvalue is a number, or when xvalue is
     a string that contains a formatted number, such as "0", "23",
     "-12.3", "4.2e3", "1.e-1", ...
     
   SEE ALSO:
 */
{
  if( is_void(xvalue) ) return 0;
  if( xvalue=="" ) return 0;
  // what's the type of xvalue ?
  tv = typeof(xvalue);
  if( anyof(tv==["char","short","int","long","float","double"]) )
    return 1;
  else if( tv=="string" ) {
    buf = strchar( xvalue );
    n = numberof(buf);
    if( buf(1)=='+' || buf(1)=='-' ) {
      if( n>2 ) {
        buf=buf(2:);
        n--;
      }
      else
        return 0;
    }
    // parse digits
    i=1;
    while( (buf(i)<='9' && buf(i)>='0') && i<n ) i++;
    nb_digits = i-1;
    if( i==n )   // no more digits at the end of the string : it's a number ! (+3, 0, 12, -32 ..)
      return 1;
    // let's now examine the end of the string
    buf = buf(i:);   // cut the first part already examined
    n -= i-1;        // update the length
    if(buf(1)=='.') {
      i=2;
      while( (buf(i)<='9' && buf(i)>='0') && i<n ) i++;
      nb_digits += i-2;
      if( i==n ) {  // no more digits at the end of the string : it's a number ! (+3, 0, 12, -32 ..)
        if( nb_digits )
          return 1;
        else
          return 0;
      } else {
        // there are some characters left ... now it should only be something like e+03 or so
        buf = buf(i:);   // cut the part already examined
        n -= i-1;        // update the length
      }
    }
    if( buf(1)=='e' ) {
      if( n>2 ) {
        buf=buf(2:);
        n--;
        if( buf(1)=='+' || buf(1)=='-' ) {
          if( n>2 ) {
            buf=buf(2:);
            n--;
          }
          else
            return 0;
        }
        // parse digits
        i=1;
        while( (buf(i)<='9' && buf(i)>='0') && i<n ) i++;
        if( i==n && nb_digits )   // no more digits at the end of the string : it's a number ! (+3, 0, 12, -32 ..)
          return 1;
        else
          return 0;
        
      } else
        return 0;
    }
    return 0;
  }
  return 0;
}





func replaceFitsKey(filename, keyword, xvalue)
/* DOCUMENT replaceFitsKey, filename, keyword, value
     Replaces the value assigned to a keyword, by a new value.
   SEE ALSO:
 */
{
  if( NO_FILE_WRITING==1 ) return;
  fic = open( filename,"r+b" );
  tmp = array(char,80);
  keywd = "";
  irec = 0;
  while ( keywd!=keyword && irec<100000 && keywd!="END" ) {
    irec++;
    _read, fic, (irec-1)*80, tmp;
    str = string(&tmp);
    keywd = mystrtrim(tmp(1:8));
  }
  if( keywd=="END" | irec==100000 ) {
    //    write,format="Keyword %s not found at record %d.\n",keyword, irec;
    close, fic;
    system,"rm "+filename+"L > /dev/null 2>&1";
    addFitsKey, filename, keyword, xvalue, noerr=1;
  } else {
    value = anything2str(xvalue);
    if( is_a_number( xvalue ) )
      tmp = putBinaryRight(value, tmp, 11, 30);   // 11 to 30 = indexes for keyword values in FITS definition
    else
      tmp = putBinaryLeft("'"+value+"'", tmp, 11, 80);   // 11 to 30 = indexes for keyword values in FITS definition
    _write, fic, (irec-1)*80, tmp;
    close, fic;
    system,"rm "+filename+"L > /dev/null 2>&1";
  }
}




func readFitsAllKey(filename)
/* DOCUMENT res = readFitsAllKey(filename)
     Reads all keywords and their value.
     The result is a string array res=array(string,N,2)
     res(,1) is the list of keywords
     res(,2) is the list of values
     
   SEE ALSO:
 */
{
  fic = open( filename,"rb" );
  tmp = array(char,80);
  keywd = "";
  tabkeywd = [keywd];
  value = "not found";
  tabvalue = [value];
  endkey = "";
  irec = 0;
  k = 0;
  while ( irec<10000 && endkey!="END" ) {
    irec++;
    _read, fic, (irec-1)*80, tmp;
    endkey = mystrtrim(tmp(1:3));   // looking for 'END' keyword
    if( tmp(9)=='=' ) {     // this is a line with keyword !
      k++;
      keywd = mystrtrim(tmp(1:8));
      if( anyof(tmp(10:30)=='\'') )
        value = getFitsString(tmp);
      else
        value = mystrtrim(tmp(10:31));
      if( k>numberof(tabvalue) ) {
        tabvalue=grow(tabvalue,tabvalue); // doubler la taille du tableau
        tabkeywd=grow(tabkeywd,tabkeywd);
      }
      tabkeywd(k) = keywd;
      tabvalue(k) = value;
    }
  }
  close, fic;
  if(k==1)
    return [keywd,value];
  else
    return [tabkeywd(1:k), tabvalue(1:k)];
}






func readFitsKey(filename, keyword)
/* DOCUMENT value = readFitsKey(filename, keyword)
     Reads the value assigned to the keyword <keyword> in the FITS file <filename>.
     When the keyword is present several times in the file, then an array of values
     will be returned.
     
   SEE ALSO:
 */
{
  fic = open( filename,"rb" );
  tmp = array(char,80);
  keywd = "";
  value = keyword + " not found";
  tabvalue = [value];
  
  irec = 0;
  k = 0;
  while ( irec<100000 && keywd!="END" ) {
    irec++;
    _read, fic, (irec-1)*80, tmp;
    str = string(&tmp); 
    keywd = mystrtrim(tmp(1:8));
    if( keywd==keyword ) {
      k++;
      if( k>numberof(tabvalue) ) tabvalue=grow(tabvalue,tabvalue); // doubler la taille du tableau
      if( keyword=="COMMENT" ) {
        // COMMENT case
        value = mystrtrim(tmp(9:80));
      } else {
        // general case
        if( anyof(tmp(10:30)=='\'') )
          value = getFitsString(tmp);
        else
          value = mystrtrim(tmp(11:30));
      }
      tabvalue(k) = value;
    }
  }
  close, fic;
  if(k<=1)
    return value;
  else
    return tabvalue(1:k);
}




func getFitsString(tmp)
{
  quotes = where(tmp=='\'');
  if( numberof(quotes)<2 )
    return  mystrtrim(tmp(11:30));
  deb = min(quotes);
  fin = max(quotes);
  if( fin<=deb )
    return  mystrtrim(tmp(11:30));
  if( fin==deb+1 )
    return "";
  return strchar( tmp(deb+1:fin-1) );

}





func mystrtrim( x )
/* DOCUMENT str = mystrtrim( charray )

   Returns a string without leading and/or trailing whitespaces, from
   an input array of char.
     
   SEE ALSO:
 */
{
  if( x==string(0) || x=="" ) return "";
  if( typeof(x)!="char" )
    return "";
  n = numberof(x);
  if( n==0 ) return "";
  if( n==1 ) x=[x];
  while( x(n)==' ' & n>1 ) {
    n--;
  }
  x = x(1:n);
  n = numberof(x);
  if( n==1 ) x=[x];
  if( n==0 ) return "";
  i=1;
  while( x(i)==' ' & i<n ) {
    i++;
  }
  return strchar(x(i:));
}



func addFitsKey(filename, keyword, xvalue, noerr=)
/* DOCUMENT addFitsKey, filename, keyword, value
     Adds a keyword assigned to a given value, in the FITS header of the file <filename>.
     The keyword is added at the end of the header, just before the END statement.

     Example
     > addFitsKey, "tmp.fits", "DATE", "30-10-1967"
     
   SEE ALSO:
 */
{
  if( NO_FILE_WRITING==1 ) return;
  if( noerr==1 ) {
    if( is_void(xvalue) ) return;
    if( is_void(filename) ) return;
    if( is_void(keyword) ) return;
  }
  // go !
  fic = open( filename,"r+b" );
  if( noerr==1 )
    if( fic==[] ) return;
  tmp = array(char,80);
  keywd = "";
  MAXREC=100000;
  
  // finding the keyword END or COMMENT.
  irec = 0;
  k = 0;
  icomment = 0;
  while ( irec<MAXREC && keywd!="END" ) {
    irec++;
    _read, fic, (irec-1)*80, tmp;
    str = string(&tmp); 
    keywd = mystrtrim(tmp(1:8));
    if( keywd==COMMENT & icomment==0 )  // des qu'un COMMENT est trouve, on repere le record
      icomment = irec;
  }
  
  if( irec==MAXREC ) {
    write,format="Unable to add keyword : the end of the FITS header has not been found after %d 80-bytes records\n",MAXREC;
    write,"Nothing done. File is not modified.";
    close, fic;
    return;
  }

  if( (irec%36)==0 ) {             // one header block (2880 bytes) has to be added
    system,"cp "+filename+" "+tmpDir+"local_addFitsKey_copy.fits";
    system_chmod,666, tmpDir+"local_addFitsKey_copy.fits";

    sfic = sizeof(fic);             // determine file size
    buf = array(char,sfic);         // preparing a big array to read the whole file
    nread = _read(fic, 0, buf);     // reads all file as a whole
    if( nread!=sfic ) {             // an error occured ...
      write,format="An error occured. Number of bytes read (%d) is different from file size (%d).\n",nred,sfic;
      write,format="Exiting before it's too late. File has not been modified.";
      close, fic;
    }
    buf = grow(buf, array(char,2880));                     // adding 2880 bytes to the buffer
    buf( irec*80+1+2880:sfic+2880 ) = buf(irec*80+1:sfic); // shifting the whole data by 2880
    buf( irec*80+1:irec*80+2880 ) = ' ';                   // filling in the new block space with ' '
    _write, fic, 0, buf;                                   // writing the new file with added header block
  }

  // writing the new keyword+value
  value = anything2str(xvalue);
  tmp = array(' ',80*2);
  tmp(9)='=';
  tmp(32)='/';
  tmp = putBinaryLeft(keyword, tmp, 1, 8);
  if( is_a_number( xvalue ) ) {
    tmp = putBinaryRight(value, tmp, 11, 30);
    // getting system date & hour
    pfic = popen("date \"+%Y-%m-%d_%Hh%Mm%Ss\"", 0);
    systemDate = rdline(pfic);
    close, pfic;
    tmp = putBinaryLeft("added keyword "+systemDate, tmp, 34, 80);
  } else {
    tmp = putBinaryLeft("'"+value+"'", tmp, 11, 80);
  }
  // writing the END
  tmp = putBinaryLeft("END", tmp, 81, 83);
  _write, fic, (irec-1)*80, tmp;                   // writing the new keyword + END
  
  close, fic;
  system,"rm "+filename+"L > /dev/null 2>&1";
}




func addFitsAllKey(filename, keywordList)
/* DOCUMENT addFitsAllKey, filenamepath, keywordList;
   
     addFitsAllKey, path, [["NBF","DTIME"],["33","0.47893"]]
     
   SEE ALSO:
 */
{
  if( NO_FILE_WRITING==1 ) return;
  n = dimsof(keywordList)(2);   // numberof keywords
  for(i=1;i<=n;i++) {
    replaceFitsKey, filename, keywordList(i,1), keywordList(i,2);
  }
}




func readFitsComment(filename)
{
  tmp = readFitsKey(filename,"COMMENT");
  // suppress empty trailing lines
  n = numberof(tmp);
  if( n<2 )
    return tmp;
  while( n>1 & tmp(n)==string(0) ) {
    n--;
  }
  return tmp(:n);
}




func addFitsComment(filename, value)
/* DOCUMENT addFitsComment, filename, value
     
   SEE ALSO:
 */
{
  if( NO_FILE_WRITING==1 ) return;
  fic = open( filename,"r+b" );
  tmp = array(char,80);
  keywd = "";
  MAXREC=100000;
  nbcom = numberof(value);   // number of comment lines to be added
  
  // finding the keyword END or COMMENT.
  irec = 0;
  k = 0;
  while ( irec<MAXREC && keywd!="END" ) {
    irec++;
    _read, fic, (irec-1)*80, tmp;
    str = string(&tmp); 
    keywd = mystrtrim(tmp(1:8));
  }
  
  if( irec==MAXREC ) {
    write,format="Unable to add keyword : the end of the FITS header has not been found after %d 80-bytes records\n",MAXREC;
    write,"Nothing done. File is not modified.";
    close, fic;
    return;
  }


  nblock = (irec-1+nbcom)/36 - (irec-1)/36;
  
  if( nblock>0 ) {             // some header block (2880 bytes) has to be added
    system,"cp "+filename+" "+tmpDir+"local_addFitsKey_copy.fits";
    system_chmod,666,tmpDir+"local_addFitsKey_copy.fits";

    sfic = sizeof(fic);             // determine file size
    buf = array(char,sfic);         // preparing a big array to read the whole file
    nread = _read(fic, 0, buf);     // reads all file as a whole
    if( nread!=sfic ) {             // an error occured ...
      write,format="An error occured. Number of bytes read (%d) is different from file size (%d).\n",nred,sfic;
      write,format="Exiting before it's too late. File has not been modified.";
      close, fic;
    }
    BLMEM = 2880*nblock;                                     // amount of space multiple of 2880 added
    buf = grow(buf, array(char,BLMEM));                      // adding k*2880 bytes to the buffer
    buf( irec*80+1+BLMEM:sfic+BLMEM ) = buf(irec*80+1:sfic); // shifting the whole data by 2880
    buf( irec*80+1:irec*80+2880 ) = ' ';                     // filling in the new block space with ' '
    _write, fic, 0, buf;                                     // writing the new file with added header block
  }

  // writing the new keyword+value
  tmp = array(' ',80);
  for(i=1;i<=nbcom;i++) {
    tmp = putBinaryLeft("COMMENT", tmp, 1, 9);
    tmp = putBinaryLeft(value(i), tmp, 9, 80);
    _write, fic, (irec-1+i-1)*80, tmp;                   // writing the new COMMENT
  }
  // writing the END
  tmp = putBinaryLeft("END", tmp, 1, 80);
  _write, fic, (irec-1+nbcom)*80, tmp;                   // writing the END keyword
  
  close, fic;
  system,"rm "+filename+"L > /dev/null 2>&1";
}




func replaceFitsComment( filename, value )
{
  if( NO_FILE_WRITING==1 ) return;
  fic = open( filename,"r+b" );
  tmp = array(char,80);
  keywd = "";
  MAXREC=100000;
  nbcom = numberof(value);   // number of comment lines to be added
  
  // finding the keyword END or COMMENT.
  irec = 0;
  k = 0;
  while ( irec<MAXREC && keywd!="END" ) {
    irec++;
    _read, fic, (irec-1)*80, tmp;
    str = string(&tmp); 
    keywd = mystrtrim(tmp(1:8));
  }
  
  if( irec==MAXREC ) {
    write,format="Unable to add keyword : the end of the FITS header has not been found after %d 80-bytes records\n",MAXREC;
    write,"Nothing done. File is not modified.";
    close, fic;
    return;
  }

  // re-parsing the header again
  irec = 0;
  k = 0;
  keywd="";
  while ( irec<MAXREC && keywd!="END" ) {
    irec++;
    _read, fic, (irec-1)*80, tmp;
    str = string(&tmp); 
    keywd = mystrtrim(tmp(1:8));
    if( keywd=="COMMENT" ) {
      k++;
      tmp = putBinaryLeft("COMMENT", tmp, 1, 80);
      if( k>numberof(value) )
        tmp = putBinaryLeft(" ", tmp, 9, 80);
      else
        tmp = putBinaryLeft(value(k), tmp, 9, 80);
      _write, fic, (irec-1)*80, tmp;                   // writing the new COMMENT
    }
  }
  close, fic;
  system,"rm "+filename+"L > /dev/null 2>&1";
 
  // writing the complement that has not been written yet
  if( k<numberof(value) )
    addFitsComment, filename, value(k+1:);
}

