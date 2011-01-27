class bedfeature(object):
    def __init__(self, chr,start,stop,
                 name=None,value=None,strand=None,
                 thickStart=None,thickStop=None,itemRGB=None,
                 blockCount=None,blockSizes=None,blockStarts=None):
        self.chr=chr
        self.start=int(start)
        self.stop=int(stop)
        self.name=name
        self.strand=strand
        try:
            self.value=float(value)
        except (TypeError,ValueError):
            self.value=value
        try:
            self.thickStart=int(thickStart)
        except TypeError:
            self.thickStart=thickStart
        try:
            self.thickStop=int(thickStop)
        except TypeError:
            self.thickStop=thickStop
        try:
            self.blockCount=int(blockCount)
        except TypeError:
            self.blockCount=blockCount
        
        self.itemRGB=itemRGB
        self.blockSizes=blockSizes
        self.blockStarts=blockStarts

    def __repr__(self):
        return 'bed feature: %s:%s-%s' % (self.chr,self.start,self.stop)
    
    def __len__(self):
        return self.stop - self.start

    def tostring(self,fields=None):
        """Prints the bed record suitable for writing to file, newline included.
        
        In the interest of speed, does not do error-checking.
        """
        
        items = [self.chr, 
                 self.start, 
                 self.stop, 
                 self.name, 
                 self.value, 
                 self.strand, 
                 self.thickStart,
                 self.thickStop, 
                 self.itemRGB,
                 self.blockCount,
                 self.blockSizes, 
                 self.blockStarts]
        printables = []
        if fields is None:
            fields = len(items)
        for item in items[0:fields]:
            if item is None:
                printables.append('')
            else:
                printables.append(str(item))
        return '\t'.join(printables).rstrip()+'\n'

