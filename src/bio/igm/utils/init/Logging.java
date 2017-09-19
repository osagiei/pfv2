/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package bio.igm.utils.init;

import java.io.IOException;
import java.util.logging.FileHandler;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.logging.SimpleFormatter;


/**
 *
 * @author osagie
 */
public class Logging {

    public static FileHandler fh;
    public static SimpleFormatter formatter;
    private String classname;

    public Logging(String path, String classname)
            throws IOException {
        fh = new FileHandler(path + "/run.log", true);
        
        formatter = new SimpleFormatter();
        this.classname = classname;
    }

    public Logger setup() {
        Logger logger = Logger.getLogger(this.classname);
        logger.setLevel(Level.ALL);
        fh.setFormatter(formatter);
        logger.addHandler(fh);
        logger.setUseParentHandlers(false);
        return logger;
    }
}
