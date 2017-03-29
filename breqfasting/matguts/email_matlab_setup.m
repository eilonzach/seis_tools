function email_matlab_setup(email_address,outgoing_smtp_server)
% email_matlab_setup(email_address,outgoing_smtp_server)
% 
% This function sets MATLAB's parameters to allow you to send emails from
% the MATLAB command line.
% 
% INPUTS:
% email_address - string email address  e.g. 'Name_Surname@brown.edu'
% outgoing_smtp_server - string outgoing server (probably 'smtp.gmail.com') 
% 
% N.B. if you don't give it any inputs, it will only set the properties to
% allow emails, not the preferences you may need to send them. This usage
% assumes you've already set your Internet prefs accordingly - try 
%  getpref('Internet') to test this. 
% 
% Z. Eilon (2015)

if nargin <2
    email_address = getpref('Internet','E_mail');
end

props = java.lang.System.getProperties;
% props.setProperty('mail.smtp.port','587');
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');

if nargin > 1
password = inputdlg('type email password here: '); password = password{1};
% Apply prefs and props
setpref('Internet','E_mail',email_address);
setpref('Internet','SMTP_Server',outgoing_smtp_server);
setpref('Internet','SMTP_Username',email_address);
setpref('Internet','SMTP_Password',password);
% Send test email
sendmail(email_address, 'Test', 'Msg from MATLAB');
end
