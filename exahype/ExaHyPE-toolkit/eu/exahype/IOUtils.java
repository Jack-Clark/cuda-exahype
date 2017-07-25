package eu.exahype;

import java.io.InputStream;
import java.util.Scanner;

public class IOUtils {
  /**
   * Converts the content of a resource which might be
   * located in a jar or in folder to a string.
   * 
   * @note The relative path must not start with a "/"!
   * 
   * <h2>ExaHyPE Toolkit</h2>
   * This trick enables us to bundle the templates with
   * the class files in a single jar.
   * 
   * <h2>Further information</h2>
   * Scanner trick: 
   * http://web.archive.org/web/20140531042945/https://weblogs.java.net/blog/pat/archive/2004/10/stupid_scanner_1.html
   * 
   * Class loader behaviour:
   * http://stackoverflow.com/questions/13269556/strange-behavior-of-class-getresource-and-classloader-getresource-in-executa
   */
  public static String convertRessourceContentToString(String relativePath) {
    InputStream in = ClassLoader.getSystemClassLoader().getResourceAsStream(relativePath);
    java.util.Scanner s = new java.util.Scanner(in);
    Scanner sNew = s.useDelimiter("\\A");
    String content = sNew.hasNext() ? sNew.next() : "";
    sNew.close();
    s.close();
    return content;
  } 
}
