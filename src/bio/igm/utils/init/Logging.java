package bio.igm.utils.init;

import java.io.File;
import java.io.IOException;
import java.util.logging.ConsoleHandler;
import java.util.logging.FileHandler;
import java.util.logging.Formatter;
import java.util.logging.Handler;
import java.util.logging.Level;
import java.util.logging.LogRecord;
import java.util.logging.Logger;

/**
 * Sets up a logger that writes to {@code <path>/run.log} and mirrors the same
 * records to stderr so that pipeline progress is visible in the shell.
 *
 * @author osagie
 */
public class Logging {

    private static final Formatter FORMATTER = new Formatter() {
        @Override
        public String format(LogRecord record) {
            StringBuilder sb = new StringBuilder(160);
            sb.append(String.format("%1$tF %1$tT", record.getMillis()))
              .append(" [").append(record.getLevel().getName()).append("] ")
              .append(shortName(record.getLoggerName())).append(" - ")
              .append(formatMessage(record)).append(System.lineSeparator());
            if (record.getThrown() != null) {
                sb.append(record.getThrown()).append(System.lineSeparator());
                for (StackTraceElement e : record.getThrown().getStackTrace()) {
                    sb.append("\tat ").append(e).append(System.lineSeparator());
                }
            }
            return sb.toString();
        }

        private String shortName(String name) {
            if (name == null) {
                return "";
            }
            int dot = name.lastIndexOf('.');
            return dot < 0 ? name : name.substring(dot + 1);
        }
    };

    private final String classname;
    private final String path;

    /**
     * @param path      directory the run log is written to; created if absent
     * @param classname logger name, normally the calling class
     * @throws IOException if the log file cannot be opened
     */
    public Logging(String path, String classname) throws IOException {
        this.path = path;
        this.classname = classname;
        File dir = new File(path);
        if (!dir.isDirectory() && !dir.mkdirs()) {
            throw new IOException("Cannot create log directory: " + dir.getAbsolutePath());
        }
    }

    /**
     * Builds the logger.  Safe to call more than once for the same class: any
     * handlers attached by a previous call are closed and replaced, so repeated
     * setup does not leak file descriptors or duplicate every line.
     *
     * @return a configured logger
     * @throws IOException if the run log cannot be opened
     */
    public Logger setup() throws IOException {
        Logger logger = Logger.getLogger(this.classname);
        for (Handler existing : logger.getHandlers()) {
            logger.removeHandler(existing);
            existing.close();
        }

        FileHandler fileHandler = new FileHandler(new File(path, "run.log").getPath(), true);
        fileHandler.setFormatter(FORMATTER);
        fileHandler.setLevel(Level.ALL);

        ConsoleHandler consoleHandler = new ConsoleHandler();
        consoleHandler.setFormatter(FORMATTER);
        consoleHandler.setLevel(Level.INFO);

        logger.setLevel(Level.ALL);
        logger.addHandler(fileHandler);
        logger.addHandler(consoleHandler);
        logger.setUseParentHandlers(false);
        return logger;
    }

    /**
     * Convenience factory: resolves {@code target} to a directory (using its
     * parent when a file is given) and returns a logger writing there.  Falls
     * back to a plain stderr logger if the run log cannot be opened, so that
     * logging never takes the pipeline down.
     *
     * @param target    working directory, or a file inside it
     * @param clazz     the calling class
     * @return a usable logger, never null
     */
    public static Logger forWorkingDir(String target, Class<?> clazz) {
        File f = new File(target);
        String dir = f.isDirectory() ? target : (f.getParent() == null ? "." : f.getParent());
        try {
            return new Logging(dir, clazz.getName()).setup();
        } catch (IOException ex) {
            Logger fallback = Logger.getLogger(clazz.getName());
            fallback.log(Level.WARNING, "Could not open run.log in " + dir + "; logging to stderr only", ex);
            return fallback;
        }
    }
}
