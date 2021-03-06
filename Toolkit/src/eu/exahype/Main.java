package eu.exahype;

public class Main {
  
  private static String CLEAN_OPT_KERNEL = "--clean-opt";
  private static String NOT_INTERACTIVE  = "--not-interactive";
  private static String INTERACTIVE      = "--interactive";
  
  public static void printHeader() {
    System.out.println("================================");
    System.out.println(" ___          _  _      ___ ___");
    System.out.println("/ __|_ ____ _| || |_  _/ _ \\ __|");
    System.out.println("| _|\\ \\ / _` | __ | || |  _/ _| ");
    System.out.println("\\___/_\\_\\__,_|_||_|\\_, |_| \\___|");
    System.out.println("                   |__/         ");
    System.out.println("================================");
    System.out.println("");
    System.out.println(" www.exahype.eu ");
    System.out.println("");
    System.out.println("================================");
    System.out.println("");
    System.out.println("The project has received funding from the European Union's ");
    System.out.println("Horizon 2020 research and innovation programme under grant ");
    System.out.println("agreement No 671698 (ExaHyPE). It is based upon the PDE ");
    System.out.println("framework Peano (www.peano-framework.org).");
    System.out.println("");
    System.out.println("");
  }

  public static void waitForInteraction(boolean interactive) {
    if (interactive) {
      System.out.println("<press Enter>");
      try {
        System.in.read();
      } catch (Exception e) {
      }
      for (int i = 0; i < 50; ++i) System.out.println();
      printHeader();
    }
  }

  /*
  // expose global variables.
  // This is ugly. Main should be a singleton instead.
  public static boolean interactive;
  public static String inputFileName;
  */

  public static void main(String[] args) {
    //
    // Usually, I write the header directly before a new algorithm phase, but
    // not for the first phase
    //    
    printHeader();

    if (args.length != 1 && args.length != 2) {
      System.err.println("ERROR: Please provide input file as argument");
      System.exit(-1);
    };

    if (args.length == 2 && args[0].compareTo(Main.NOT_INTERACTIVE) != 0
        && args[0].compareTo(Main.INTERACTIVE) != 0) {
      System.err.println(
          "ERROR: First optional argument has to be "+Main.NOT_INTERACTIVE+" or "+Main.INTERACTIVE+". Received \""
          + args[0] + "\"");
      System.exit(-2);
    };

    boolean interactive = args.length == 1 || args[0].compareTo(Main.INTERACTIVE) == 0;

    if (args.length == 1) {
      System.out.println(
          "INFO: You might want to add "+Main.NOT_INTERACTIVE+" or "+Main.INTERACTIVE+" as first command ");
      System.out.println("      line argument to control whether script runs interactively");
    }

    //
    // Parse file
    //
    String inputFileName = args.length == 2 ? args[1] : args[0];

    eu.exahype.parser.Parser parser = null;
    eu.exahype.node.Start document = null;
    try {
      System.out.print("read input file " + inputFileName + " .");
      parser = new eu.exahype.parser.Parser(new eu.exahype.lexer.Lexer(
          new java.io.PushbackReader(new java.io.FileReader(inputFileName),5000))); //TODO Dominic fix the PushbackReader buffer size
      document = parser.parse();
      System.out.println(".. ok");
      System.out.println("\n\n\n\n");
      System.out.println("Start to interpret script ... ");
      waitForInteraction(interactive);
    } catch (Exception e) {
      System.out.println(".. failed ");
      System.out.println("\n\n\n\n");
      System.out.println("ERROR: " + e.toString());
      System.exit(-3);
    }

    DirectoryAndPathChecker directoryAndPathChecker = null;

    //
    // Check directories and pathes
    //
    try {
      directoryAndPathChecker = new DirectoryAndPathChecker();

      document.apply(directoryAndPathChecker);

      System.out.println("\n\n\n\n");
      if (!directoryAndPathChecker.valid) {
        System.err.println("ERROR: Some directories did not exist and/or could not be created");
        System.err.println("ExaHyPE script failed ");
        System.exit(-4);
      }
      System.out.println("validated and configured pathes ... ok");
      waitForInteraction(interactive);
    } catch (Exception e) {
      System.out.println("ERROR: " + e.toString());
      System.err.println("ExaHyPE script failed ");
      System.exit(-5);
    }
    
    // Create the solvers
    try {
      CreateSolverClasses createSolverClasses = new CreateSolverClasses(directoryAndPathChecker);

      document.apply(createSolverClasses);

      System.out.println("\n\n\n\n");
      if (!createSolverClasses.valid) {
        System.err.println("ERROR: Could not create application's solver classes");
        System.err.println("ExaHyPE script failed ");
        System.exit(-6);
      }
      System.out.println("generate application-specific solver classes ... ok");
      waitForInteraction(interactive);
    } catch (Exception e) {
      System.out.println("ERROR: " + e.toString());
      // @tood remove later again
      e.printStackTrace();
      System.err.println("ExaHyPE script failed ");
      System.exit(-7);
    }

    // Create the plotters
    try {
      CreatePlotterClasses createPlotterClasses = new CreatePlotterClasses(directoryAndPathChecker);

      document.apply(createPlotterClasses);

      System.out.println("\n\n\n\n");
      if (!createPlotterClasses.valid) {
        System.err.println("ERROR: Could not create application's plotter classes");
        System.err.println("ExaHyPE script failed ");
        System.exit(-8);
      }
      System.out.println("generate application-specific plotter classes ... ok");
      waitForInteraction(interactive);
    } catch (Exception e) {
      System.out.println("ERROR: " + e.toString());
      System.err.println("ExaHyPE script failed ");
      e.printStackTrace();
      System.exit(-9);
    }

    // Create the kernel calls
    try {
      GenerateSolverRegistration generateKernelCalls =
          new GenerateSolverRegistration(directoryAndPathChecker, inputFileName);

      document.apply(generateKernelCalls);

      System.out.println("\n\n\n\n");
      if (!generateKernelCalls.valid) {
        System.err.println("ERROR: Could not create ExaHyPE's kernel calls");
        System.err.println("ExaHyPE script failed ");
        System.exit(-10);
      }
      System.out.println("generate computational kernel calls ... ok");
      waitForInteraction(interactive);
    } catch (Exception e) {
      System.out.println("ERROR: " + e.toString());
      System.err.println("ExaHyPE script failed ");
      e.printStackTrace();
      System.exit(-11);
    }

    //
    // Setup build environment, i.e. makefiles
    //
    try {
      SetupBuildEnvironment setupBuildEnvironment =
          new SetupBuildEnvironment(directoryAndPathChecker);

      document.apply(setupBuildEnvironment);

      System.out.println("\n\n\n\n");
      if (!setupBuildEnvironment.valid) {
        System.err.println("ERROR: Could not create ExaHyPE's build environment");
        System.err.println("ExaHyPE script failed ");
        System.exit(-12);
      }
      System.out.println("setup build environment ... ok");
      waitForInteraction(interactive);
    } catch (Exception e) {
      System.out.println("ERROR: " + e.toString());
      System.err.println("ExaHyPE script failed ");
      e.printStackTrace();
      System.exit(-13);
    }
  }
}
