// These tools are for semi-automated manual length measurement of cytonemes length in a time lapse or static image
// The user may use these tools to measure cytonemes/filopodia length manually basically by drawing a freehand line 
// over the same cytoneme over time or different cytonemes in the same static image
// Run these tools with FIJI
//
// Authors: Guilherme Oliveira Barbosa and Thomas B. Kornberg
// Filliations: University of California San Francisco - USA
// Support contact: gbarbosa.bio@gmail.com
//
//2021

//Starting macro bottom and dialog box
macro "Start CytoID Tool - N66C000D00D01D02D03D04D05D06D07D08D09D0aD0bD0cD0dD0eD10D11D12D13D14D15D16D17D18D19D1aD1bD1cD1dD1eD20D21D22D23D24D25D26D27D28D29D2aD2bD2cD2dD2eD30D31D32D33D34D35D36D37D38D39D3aD3bD3cD3dD3eD40D41D42D43D44D45D46D47D48D49D4aD4bD4cD4dD4eD50D51D52D5aD5bD5cD5dD5eD60D61D62D66D6aD6bD6cD6dD6eD70D71D72D74D75D76D77D78D7aD7bD7cD7dD7eD80D81D82D85D86D87D8aD8bD8cD8dD8eD90D91D92D95D96D97D9aD9bD9cD9dD9eDa0Da1Da2Da3Da4Da5Da6Da7Da8Da9DaaDabDacDadDaeDb0Db1Db2Db3Db4Db8Db9DbaDbcDbdDbeDc0Dc1Dc2Dc3Dc4DccDcdDceDd0Dd1Dd2Dd3Dd4Dd5Dd6DdcDddDdeDe0De1De2De3De4DebDecDedDeeC000D53D59D65D67DeaCfffD54D55D56D57D58D63D64D68D69D73D79D83D84D88D89D94D98Db5DbbDc5Dc6Dc7Dc8DcaDcbDd8Dd9DdaDdbDe6De7De8C222D99Dd7C111Db7C666C555Dc9De9C222C777Db6C333D93De5Bf0C000D00D01D02D03D04D08D09D0aD0bD0cD0dD0eD10D11D12D1aD1bD1cD1dD1eD20D21D22D23D24D26D27D28D2aD2bD2cD2dD2eD30D31D32D33D34D35D36D37D38D39D3aD3bD3cD3dD3eD40D41D42D43D44D4aD4bD4cD4dD4eD50D51D52D53D54D57D5aD5bD5cD5dD5eD60D61D62D63D64D66D67D68D6aD6bD6cD6dD6eD70D71D72D73D74D77D7aD7bD7cD7dD7eD80D81D82D83D84D8aD8bD8cD8dD8eD90D91D92D93D94D95D96D97D98D99D9aD9bD9cD9dD9eDa0Da1Da2Da3Da4Da5Da6Da7Da8Da9DaaDabDacDadDaeC000CfffD05D13D14D15D16D17D18D19D25D29D46D47D48D55D59D65D69D75D79D86D87D88C222D49C111D07C666D56C555D58D76D78C222D45D85D89C777D06C333B0fC000D01D02D03D04D05D07D08D09D0aD13D14D17D18D19D1aD27D28D29D2aD30D31D32D33D34D35D36D37D38D39D3aD40D41D42D43D44D45D46D47D48D49D4aD50D51D52D53D54D55D56D57D58D59D5aD60D61D62D63D64D65D66D67D68D69D6aD70D71D72D73D74D75D76D77D78D79D7aD80D81D82D83D84D85D86D87D88D89D8aD90D91D92D93D94D95D96D97D98D99D9aDa0Da1Da2Da3Da4Da5Da6Da7Da8Da9DaaC000D12CfffD00D06D10D11D15D16D21D22D23D24D25C222C111D20D26C666C555C222C777C333Nf0C000D00D01D02D03D04D05D06D07D08D09D0aD10D11D12D13D14D15D16D17D18D19D1aD20D21D22D23D24D25D26D27D28D29D2aD30D31D32D33D34D35D36D37D38D39D3aD40D41D42D43D44D45D46D47D48D49D4aD50D51D52D53D54D55D56D57D58D59D5aD60D61D62D63D64D65D66D67D68D69D6aD70D71D72D73D74D75D76D77D78D79D7aD80D81D82D83D84D85D86D87D88D89D8aD90D91D92D93D94D95D96D97D98D99D9aDa7Da8Da9DaaDb0Db1Db2Db3Db4Db5Db6Db7Db8Db9DbaDc0Dc1Dc2Dc3Dc4Dc5Dc6Dc7Dc8Dc9DcaDd7Dd8Dd9DdaDe1De2De3De4De5De7De8De9DeaC000CfffDa0Da1Da2Da3Da4Da5Da6Dd0Dd1Dd2Dd3Dd4Dd5Dd6De0De6C222C111C666C555C222C777C333"{
      requires("1.39r");
getDimensions(width, height, channels, slices, frames);
getPixelSize(unit, pixelWidth, pixelHeight);
Dialog.create("Welcome to Cytoneme/Filopodia analysis tool on FIJI/ImageJ");
Dialog.addMessage("Before you start your cytoneme or filopodia measurements, we need some information here. \n 1-Provide the name o the experimental group in the first box. \n 2- Provide the sample number for this image. \n 3- Calibrate the time adding how many min OR sec per frame. \n 4- Add the space that value each pixel represents and indicate the unit in the folowing box");
TIME= 1;
PIXEL = pixelWidth;
unity = newArray("pixel","cm","mm","um","nm");
Dialog.addString("Group name", "Type name");
Dialog.addNumber("Sample number", 0);
Dialog.addNumber("Time calibration", TIME);
Dialog.addNumber("Pixel calibration", PIXEL);
Dialog.addChoice("Pixel unit", unity,"pixel");
Dialog.show();
Group = Dialog.getString();
Samplen = Dialog.getNumber();
TIME=Dialog.getNumber();
PIXEL=Dialog.getNumber();
unity = Dialog.getChoice();
id=d2s(random()*100000,0);

//Setting time time frame interval
Stack.setFrameInterval(TIME) 

//Ressting pixel scale
run("Set Scale...", "distance=1 known="+PIXEL+" pixel=1 unit="+unity+"");

//Get file location information
	  path = getDirectory("image");
	  pathsave=path+Group+"_Sample_"+Samplen;

//Create CytonemeID and dataframe windows
	  title1 = "CytonemeID";
      title2 = "["+title1+"]";
      title3 = getTitle()+"_data";
      title4 = "["+title3+"]";
      run("Text Window...", "name="+title2+" width=35 height=2 menu");
      run("Text Window...", "name="+title4+" width=72 height=15 menu");
      print(title4, "ID"+ "  \t" +"Xi"+ "  \t" +"Yi"+ "  \t" +"Xf"+ "  \t" +"Yf"+ "  \t" +"Time"+ "  \t" +"Length"+"\n" );
      	print("[CytonemeID]", "\\Update:");
  		id=d2s(random()*100000,0);
  		print("[CytonemeID]", id);

 //Stop this macro for the user to make measurements 
      waitForUser("Measuring", "If you chose the Cytoneme Dynamics tool, \n find your first cytoneme and start measuring. \n Once you finish, push the key c and measure the next one. \n When satisfied push ok! \n If you chose the Cytoneme Static tool, \n just measure all the cytonemes in the image and push OK. \n Good measurements for you!");;
 
//Close all windows and save data frame in the same folder the image is stored
      wait(100);
      selectWindow("CytonemeID");
      run("Close");
      selectWindow(title3);
      run("Text...", "save=&pathsave");
     
}

//Random cytonemeID generator activated when key "c" is pressed in the keyboard
	macro "CytonemeID [c]"{
		if(isOpen("CytonemeID")){
  	print("[CytonemeID]", "\\Update:");
  	id=d2s(random()*100000,0);
  	print("[CytonemeID]", id);
  	}else{
  		title1 = "CytonemeID";
      	title2 = "["+title1+"]";
      	run("Text Window...", "name="+title2+" width=35 height=2 menu");
      	print("[CytonemeID]", "\\Update:");
  		id=d2s(random()*100000,0);
  		print("[CytonemeID]", id);
  	}}

// Dynamics analysis - botton activate analysis which is basically drawing a line over the same cytoneme over time
// Repeat the same procedure for as many cytoneme as desired
  	      macro "Cytoneme Dynamic Tool - N66C000D3dD45D47D4aD53D63D6cD71D75D76D7cD81D8dD8eD91D95D98D99D9dD9eDa1Da2Da5Da8DabDadDaeDb0Db1Db2Db4Db7DbbDbdDbeDc0Dc1Dc4Dc7Dc9DcbDceDd0Dd1Dd3Dd4Dd6Dd7Dd9DdbDdcDdeDe0De1De3De6De9DebDecDeeC002C000D74D84Cf96D0eD6eC618D55Dc8Dd8C001D7dCfffC316C001D16D73D86Da0De8CffbD00D01D02D03D04D05D06D07D08D09D0aD0bD0cD0dD14D1aD1bD1cD1dD1eD20D21D22D23D24D25D26D27D28D29D2aD2bD2cD2dD2eD30D31D32D33D34D35D36D37D38D39D3bD41D4cD51D58D66D67D68D77D78D7bD87D8aD8bDa6Db6Dd5Cc37D4dD4eD5aD89C001CfffC114DeaC000Cfb8D40De4C927DbaC517D17D80Cf65D3aD60D6dC102D83Cfa6D59C728D3eD70De2C417DbcDe7Ce46D9cDc6C215D10D93DdaCfd9Ca37D44C617D6bDddCf75D79DacC111Cf96D4bD52D69C718D11D43C407Da7Db3Dc2Cc37D88C215D82Cfc8Ca37D50D56D5cC517D46D5bD62D72DcaCf75D49D5dC103D85Cfa7D7aC828Dd2Cf65D42C216DcdCffbDc5Cb37D54Cf85D65C002D90D92C306D94C214D8cCfc8D5eD6aC928D64D96D9aC103Da9Db5DccDedCfa7De5Ce55D19D48DaaCfeaD13D57Cb37C112C728Db8Cc47Cfd9D12Cf75D7eCfb7C828D15D3cCf86D61C002Db9C306Da3Dc3C114D9bCe56CfdaD97Cf86C307Cd46D18Cfd9Cfb7C316Da4CfeaBf0C000D00D03D05D06D08D09D0bD0cD0dD10D12D15D17D18D19D1bD1cD20D22D24D27D28D29D2bD2cD2eD32D33D34D36D37D38D39D3bD3cD3dD3eD40D41D42D43D45D46D47D48D49D4bD4cD4dD4eD50D51D52D55D56D57D58D59D5bD5cD5dD5eD60D61D62D63D64D65D66D67D68D69D6aD6bD6cD6dD6eD70D71D72D73D74D75D76D77D78D79D7aD7bD7cD7dD7eD80D81D82D83D84D85D86D87D88D89D8aD8bD8cD8dD8eD90D91D92D93D94D95D96D97D98D99D9aD9bD9cD9dD9eDa0Da1Da2Da3Da4Da5Da6Da7Da8Da9DaaDabDacDadDaeC002C000D2dD30Cf96C618D1eD21C001D02D53CfffC316C001D25D5aCffbD04Cc37D13C001D26D54CfffC114C000Cfb8C927C517Cf65C102Cfa6C728D1dC417Ce46C215Cfd9Ca37C617Cf75C111Cf96C718C407D1aCc37C215D11D14D35D3aD44Cfc8Ca37D0eC517D0aCf75C103D31Cfa7C828Cf65C216D4aCffbCb37Cf85C002D01C306C214D2aCfc8C928C103Cfa7Ce55CfeaCb37C112C728Cc47Cfd9Cf75Cfb7C828D23Cf86C002C306C114Ce56CfdaCf86C307D16Cd46Cfd9Cfb7C316D07CfeaB0fC000D01D04D05D08D09D0aD10D11D14D15D17D18D19D1aD20D21D22D24D27D28D29D2aD30D31D39D3aD40D45D4aD50D53D54D55D56D57D5aD60D62D63D64D65D66D67D68D6aD70D73D74D75D76D77D7aD80D84D85D86D8aD90D91D94D95D99D9aDa0Da1Da2Da3Da4Da5Da6Da7Da8Da9DaaC002D03C000D02D12Cf96C618C001CfffD61C316D07C001CffbCc37D25C001CfffD32D33D34D35D36D37D38D41D42D43D47D48D49D51D52D58D59D69D71D78D79D81D82D87D88D89D92D93D97D98C114C000D44D46D96Cfb8C927C517Cf65C102D13Cfa6C728C417Ce46C215Cfd9Ca37C617D06Cf75C111D72Cf96C718C407D00Cc37C215Cfc8Ca37C517D26Cf75C103Cfa7C828Cf65C216CffbCb37Cf85C002C306C214Cfc8C928C103Cfa7Ce55CfeaCb37C112D83C728Cc47D16Cfd9Cf75Cfb7C828Cf86C002C306C114D23Ce56CfdaCf86C307Cd46Cfd9Cfb7C316CfeaNf0C000D23D32D37D44D45D63D64D74D75D76D77D79D86D87D92Da0Da2Da3Da4DaaDb0Db2Db3Db4Db6Db7DbaDc0Dc2Dc3Dc6Dc9DcaDd3Dd6Dd9DdaDe3De5De8De9DeaC002D30C000D89Cf96D50Dd1C618D59D5aD65C001CfffC316Db5C001D39Dc7Dd5De4CffbD00D01D02D03D04D05D06D07D08D09D0aD10D11D12D13D14D15D18D19D20D21D24D51D61D71D7aD81D8aD91Da1Da8Db1Db8Dc1De7Cc37C001D90Dc4Dd4CfffC114C000Cfb8C927D57C517D56D66D93Cf65D83D99C102D33Cfa6D22C728D26D6aC417D17D35D38Ce46D48C215Db9Dd0Dd2Cfd9D52Ca37D60D85C617D47Cf75D4aD78D94C111Cf96D1aD54D62C718D49C407D36D69D82Dd8Cc37D80C215Da5Cfc8D84Ca37D67C517D2aD43D46Cf75C103D29De6Cfa7D73C828Cf65D27C216De1CffbD16Cb37D3aD70D96Cf85C002D97D9aDa6De2C306C214Cfc8D72D88C928C103Dc5Cfa7D68Dc8Ce55D31CfeaD41D53Cb37D58C112C728D42Cc47Cfd9Cf75D55Cfb7Dd7C828D98De0Cf86D34C002C306C114Ce56D40CfdaCf86D25C307Cd46Cfd9Da7Cfb7D95Da9C316CfeaD28"{
   
//setting measurements to zero      
      TIME = Stack.getFrameInterval();
      getCursorLoc(x, y, z, flags);
      xarray=newArray();
      yarray=newArray();
      setOption("disablePopupMenu", true);

//drawing a freehand line and saving it spacial information
      while (flags!=0) { 
          getCursorLoc(x, y, z, flags); 
          xarray=Array.concat(xarray,x);
          yarray=Array.concat(yarray,y);
          xi=xarray[0];
          yi=yarray[0];
          xf=xarray[lengthOf(xarray)-1];
          yf=yarray[lengthOf(yarray)-1];
          makeSelection("freeline", xarray, yarray);
          wait(1); 
      }
      
//get the data and print in to dataframe window and move to next slice
      selectWindow("CytonemeID");
      id = getInfo("window.contents");
      data0=xi;
      data1=yi;
      data2=xf;
      data3=yf;
      data8= (getValue("Slice")-1)*TIME;
      List.setMeasurements;
      data14= getValue("Length");
      title3 = getTitle()+"_data";
      title4 = "["+title3+"]";
      print(title4,id+ "  \t" +data0+ "  \t" +data1+ "  \t" +data2+ "  \t" +data3+ "  \t" +data8+ "  \t" +data14+"\n");
      roiManager("Add");
      roiManager("Show None");
      run("Select None");
      run("Next Slice [>]");
      selectImage(getTitle());
      }

// Static analysis - botton activate analysis which is basically drawing a line over different cytonemes in this image
// Repeat the same procedure for as many cytoneme as desired
        macro "Cytoneme Static Tool - N66C000D3dD45D47D4aD53D63D6cD71D75D76D7cD81D8dD8eD91D95D98D99D9dD9eDa1Da2Da5Da8DabDadDaeDb0Db1Db2Db4Db7DbbDbdDbeDc0Dc1Dc4Dc7Dc9DcbDceDd0Dd1Dd3Dd4Dd6Dd7Dd9DdbDdcDdeDe0De1De3De6De9DebDecDeeC003C000Cf7eD59C60cD6bDddC002D7dCfffD65C209D10D93DcdDdaC001D16D73D86Da0De8CfefD00D01D02D03D04D05D06D07D08D09D0aD0bD0cD0dD14D1aD1bD1cD1dD1eD20D21D22D23D24D25D26D27D28D29D2aD2bD2cD2dD2eD30D31D32D33D34D35D36D37D38D39D3bD41D51D58D66D67D68D77D78D79D7bD87D8aD8bDa6Db6Dd5Cf08D40C002CfffC106D85C001D74D84CfafD0eD6eC80fDd2C002C40aDbcDe7Cc0fD4dD4eD5aD89C105D83Cf9fD57C70dD11D43C309D94Cb0fC208D82CfdfD61D7eC90fDbaC50bD46D5bD62DcaCf2fD18C004D90D92Cf8fD42C60dD55C309Da4Cf2bC107D9bCfbfC90fD64D96D9aC50bD17D80Cc0fD88C105Da9Db5DccDedCfafD13D6dC70eDe2C40aDb3Dc2Cf4cDe5C208CfefD4cDacCa0fD44C50cD72Cf6fC003Db9Cf7fD97C60dDc8Dd8Cf1aD5eC107DeaCfbfD3aD60C80fD15C70eDb8C409Da7Cb0fD54Ca0fD50D56D5cCf4fC004Cf9fD4bD52C309Da3Dc3Cf4bD7aC207D8cCfcfD5dDc5C70eD3eD70Cf5eC208Cf7fD19D48DaaCf8fD69Cf19C80fD3cCf2fD9cDc6Cf3bCfcfD49Cf4dD12Cf8fCf1aD6aCf7eCf19De4Bf0C000D00D03D05D06D08D09D0bD0cD0dD10D12D15D17D18D19D1bD1cD20D22D24D27D28D29D2bD2cD2eD32D33D34D36D37D38D39D3bD3cD3dD3eD40D41D42D43D45D46D47D48D49D4bD4cD4dD4eD50D51D52D55D56D57D58D59D5bD5cD5dD5eD60D61D62D63D64D65D66D67D68D69D6aD6bD6cD6dD6eD70D71D72D73D74D75D76D77D78D79D7aD7bD7cD7dD7eD80D81D82D83D84D85D86D87D88D89D8aD8bD8cD8dD8eD90D91D92D93D94D95D96D97D98D99D9aD9bD9cD9dD9eDa0Da1Da2Da3Da4Da5Da6Da7Da8Da9DaaDabDacDadDaeC003D54C000D30Cf7eC60cC002CfffC209D4aC001D25D5aCfefD04Cf08C002D26CfffC106D31C001D2dCfafC80fC002D02D53C40aCc0fD13C105Cf9fC70dC309D16Cb0fC208D11D14D35D44CfdfC90fC50bCf2fC004Cf8fC60dC309D07Cf2bC107CfbfC90fC50bCc0fC105CfafC70eC40aD1aCf4cC208D3aCfefCa0fC50cD0aCf6fC003Cf7fC60dD1eD21Cf1aC107CfbfC80fD23C70eC409Cb0fCa0fD0eCf4fC004D01Cf9fC309Cf4bC207D2aCfcfC70eD1dCf5eC208Cf7fCf8fCf19C80fCf2fCf3bCfcfCf4dCf8fCf1aCf7eCf19B0fC000D01D04D05D08D09D0aD10D11D14D15D17D18D19D1aD20D21D22D24D27D28D29D2aD30D31D39D3aD40D45D4aD50D53D54D55D56D57D5aD60D63D64D65D66D67D68D6aD70D73D74D75D76D77D7aD80D84D85D86D8aD90D91D94D95D99D9aDa0Da1Da2Da3Da4Da5Da6Da7Da8Da9DaaC003D03D46D96C000D12D62Cf7eC60cD06C002CfffC209C001CfefCf08C002D44CfffD32D33D34D35D36D37D38D41D42D43D47D48D49D51D52D58D59D69D71D78D79D81D82D87D88D89D92D93D97D98C106C001D02CfafC80fC002C40aCc0fD25C105D13D72Cf9fC70dC309Cb0fC208CfdfD61C90fC50bD26Cf2fC004Cf8fC60dC309D07Cf2bC107D23CfbfC90fC50bCc0fD16C105CfafC70eC40aD00Cf4cC208CfefCa0fC50cCf6fC003Cf7fC60dCf1aC107D83CfbfC80fC70eC409Cb0fCa0fCf4fC004Cf9fC309Cf4bC207CfcfC70eCf5eC208Cf7fCf8fCf19C80fCf2fCf3bCfcfCf4dCf8fCf1aCf7eCf19Nf0C000D23D32D37D44D45D63D64D74D75D76D77D79D86D87D92Da0Da2Da3Da4DaaDb0Db2Db3Db4Db6Db7DbaDc0Dc2Dc3Dc6Dc9DcaDd3Dd6Dd9DdaDe3De5De8De9DeaC003D30C000Cf7eC60cD47C002Dc7De4CfffD4aD78C209C001D39Dd5CfefD00D01D02D03D04D05D06D07D08D09D0aD10D11D12D13D14D15D18D19D20D21D24D25D51D61D71D7aD81D8aD91D94Da1Da8Db1Db8Dc1De7Cf08C002D90Dc4Dd4CfffC106D29De6C001D89CfafD53C80fC002C40aD17D35D38Cc0fC105D33Dc5Cf9fD41D54D62D99C70dD49C309Cb0fD58C208Da5CfdfD16D34D55C90fD57C50bD2aD43Cf2fC004D9aDe2Cf8fC60dD5aC309Db5Cf2bD84C107CfbfD50C90fC50bD56D66D93Cc0fD80C105CfafD83C70eD26D6aC40aD36D69Dd8Cf4cC208CfefCa0fD60D85C50cD46Cf6fD31D52C003Cf7fC60dD59D65Cf1aC107CfbfDd1C80fD98De0C70eD42C409D82Cb0fD3aD70D96Ca0fD67Cf4fD40C004D97Da6Cf9fD1aD28C309Cf4bD68D73Dc8C207CfcfC70eCf5eDa7C208Db9Dd0Dd2De1Cf7fCf8fCf19D95Da9C80fCf2fD48Cf3bDd7CfcfCf4dCf8fD27Cf1aD72D88Cf7eD22Cf19"{
      
//setting measurements to zero   
      getCursorLoc(x, y, z, flags);
      xarray=newArray();
      yarray=newArray();
      setOption("disablePopupMenu", true);

//drawing a freehand line and saving it spacial information
      while (flags!=0) { 
          getCursorLoc(x, y, z, flags); 
          xarray=Array.concat(xarray,x);
          yarray=Array.concat(yarray,y);
          xi=xarray[0];
          yi=yarray[0];
          xf=xarray[lengthOf(xarray)-1];
          yf=yarray[lengthOf(yarray)-1];
          makeSelection("freeline", xarray, yarray);
          wait(1); 
      }

//get the data and print in to dataframe window and move to next slice
      id=getTitle();
      data0=xi;
      data1=yi;
      data2=xf;
      data3=yf;
      List.setMeasurements;
      data14= getValue("Length");
      title3 = getTitle()+"_data";
      title4 = "["+title3+"]";
      print(title4,id+ "  \t" +data0+ "  \t" +data1+ "  \t" +data2+ "  \t" +data3+ "  \t" +data14+"\n");
      roiManager("Add");
      roiManager("Show None");
      run("Select None");
      selectImage(getTitle());
      }

