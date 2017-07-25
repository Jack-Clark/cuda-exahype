package eu.exahype.solvers;

import java.io.IOException;

public class IncludeOnceHelper {
	/**
	 * A class with chainable methods to represent C++ preprocessor ifndefs.
	 * This class is part of the initiative for the toolkit to generate
	 * high quality, resuable C++ code instead of spaghetti code.
	 **/
	String headername;
	java.io.BufferedWriter writer;
	
	IncludeOnceHelper(java.io.BufferedWriter _writer, String _headername) {
		writer = _writer;
		headername = _headername;
	}
	
	IncludeOnceHelper open()  throws IOException {
		writer.write("#ifndef __"+headername+"__\n");
		writer.write("#define __"+headername+"__\n");
		return this;
	}
	
	IncludeOnceHelper close()  throws IOException {
		writer.write("#endif /* __"+headername+"__ */\n");
		return this;
	}
}