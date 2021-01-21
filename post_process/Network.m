classdef Network < PostP
    %Network force chain networks
    
    properties
        Z       % Coordination number
    end
    
    methods
        
        function obj = Network(dirPath)
            obj = obj@PostP(dirPath);
        end
        
        function obj = getZ(obj)
            % Get the time averaged Z
            n = PostP.nFiles(obj.dirPath);
            nSample = 50;
            fprintf('Sample size for coordination number: %d\n', nSample);
            i = 1;
            for filenum = n:-1:n-nSample
                fp = strcat(obj.dirPath, '/contact/contact_', string(filenum));
                g = obj.readGraph(fp);
                ZZ(i) = mean(g.degree);
                i = i+1;
            end
%             plot(ZZ(ZZ>0), 'x')
            obj.Z = mean(ZZ(ZZ>0));
        end
        
        function plotdegreeDist(obj, filenum)
            fp = strcat(obj.dirPath, '/contact/contact_', string(filenum));
            g = obj.readGraph(fp);
            histogram(g.degree(), 'BinEdges', [0 1 2 3 4 5 6 7 8]-.5, 'Normalization', 'pdf');
            ylabel('pdf')
            xlabel('Node Degree')
        end
        
    end
    
    methods (Static)
        function g = readGraph(fp)
            % Read in contacts from contact file and make graph        
            contact = Contact();
            contact = readContact(contact, fp);
            g = graph(contact.ID_A(:)+1, contact.ID_B(:)+1);
        end
    end
end

